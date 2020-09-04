#include <chrono>
#include <mutex>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/async_input_buffer.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>

#include <minimiser_model.hpp>

class sync_out
{
public:
    sync_out() = default;
    sync_out(sync_out const &) = default;
    sync_out & operator=(sync_out const &) = default;
    sync_out(sync_out &&) = default;
    sync_out & operator=(sync_out &&) = default;
    ~sync_out() = default;

    sync_out(std::filesystem::path const & path) : file(std::ofstream{path}) {}

    void write(std::string const & data)
    {
        std::lock_guard<std::mutex> lock(write_mutex);
        file << data;
    }

private:
    std::ofstream file;
    std::mutex write_mutex;
};

std::vector<size_t> compute_simple_model(cmd_arguments const & arguments)
{
    std::vector<size_t> precomp_thresholds;

    if (arguments.threshold == 0.0 && !do_cerealisation_in(precomp_thresholds, arguments))
    {
        precomp_thresholds = precompute_threshold(arguments.pattern_size,
                                                  arguments.window_size,
                                                  arguments.kmer_size,
                                                  arguments.errors,
                                                  arguments.tau);

        do_cerealisation_out(precomp_thresholds, arguments);
    }

    return precomp_thresholds;
}

template <bool compressed>
struct ibf_builder
{
    cmd_arguments const * const arguments;

    auto construct()
    {
        assert(arguments != nullptr);

        seqan3::ibf_config cfg{seqan3::bin_count{64u},
                               seqan3::bin_size{1024u},
                               seqan3::hash_function_count{2u},
                               1u};
        auto view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments->kmer_size},
                                                  seqan3::window_size{arguments->window_size},
                                                  seqan3::seed{adjust_seed(arguments->kmer_size)});
        return seqan3::technical_binning_directory{std::vector<std::vector<seqan3::dna4>>{},
                                                   std::move(view),
                                                   cfg,
                                                   true};
    }

    auto ibf()
        requires !compressed
    {
        return construct();
    }

    auto ibf()
        requires compressed
    {
        auto tmp = construct();

        return seqan3::technical_binning_directory<seqan3::data_layout::compressed,
                                                   typename decltype(tmp)::hash_adaptor_t,
                                                   seqan3::dna4>{std::move(tmp)};
    }
};

template <typename t>
void load_ibf(t & tbd, cmd_arguments const & arguments, size_t const part, double & ibf_io_time)
{
    std::filesystem::path ibf_file{arguments.ibf_file};
    ibf_file += "_" + std::to_string(part);
    std::ifstream is{ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    auto start = std::chrono::high_resolution_clock::now();
    iarchive(tbd);
    auto end = std::chrono::high_resolution_clock::now();
    ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

template <typename t>
inline void do_parallel(t && worker, size_t const num_records, size_t const threads, double & compute_time)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<decltype(std::async(std::launch::async, worker, size_t{}, size_t{}))> tasks;
    size_t const records_per_thread = num_records / threads;

    for (size_t i = 0; i < threads; ++i)
    {
        size_t const start = records_per_thread * i;
        size_t const end = i == (threads-1) ? num_records: records_per_thread * (i+1);
        tasks.emplace_back(std::async(std::launch::async, worker, start, end));
    }

    for (auto && task : tasks)
        task.wait();

    auto end = std::chrono::high_resolution_clock::now();
    compute_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

template <bool compressed>
void run_program_multiple(cmd_arguments const & arguments)
{
    using std::get;

    ibf_builder<compressed> builder{&arguments};
    auto tbd = builder.ibf();

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    size_t const kmers_per_window = arguments.window_size - arguments.kmer_size + 1;
    size_t const kmers_per_pattern = arguments.pattern_size - arguments.kmer_size + 1;
    size_t const min_number_of_minimisers = kmers_per_window == 1 ? kmers_per_pattern :
                                                std::ceil(kmers_per_pattern / static_cast<double>(kmers_per_window));
    size_t const max_number_of_minimisers = arguments.pattern_size - arguments.window_size + 1;
    std::vector<size_t> const precomp_thresholds = compute_simple_model(arguments);

    auto cereal_worker = [&] ()
    {
        load_ibf(tbd, arguments, 0, ibf_io_time);
    };

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        auto cereal_handle = std::async(std::launch::async, cereal_worker);

        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        std::vector<seqan3::counting_vector<uint16_t>> counts(records.size(),
                                                              seqan3::counting_vector<uint16_t>(tbd.bin_count(), 0));

        auto count_task = [&](size_t const start, size_t const end)
        {
            auto counter = tbd.template counting_agent<uint16_t>();
            size_t counter_id = start;

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
            {
                auto & result = counter.count_query(seq);
                counts[counter_id++] += result;
            }
        };

        cereal_handle.wait();
        do_parallel(count_task, records.size(), arguments.threads, compute_time);

        for (size_t const part : std::views::iota(1u, static_cast<unsigned int>(arguments.parts - 1)))
        {
            load_ibf(tbd, arguments, part, ibf_io_time);
            do_parallel(count_task, records.size(), arguments.threads, compute_time);
        }

        load_ibf(tbd, arguments, arguments.parts - 1, ibf_io_time);
        sync_out synced_out{arguments.out_file};

        auto output_task = [&](size_t const start, size_t const end)
        {
            auto counter = tbd.template counting_agent<uint16_t>();
            size_t counter_id = start;
            std::string result_string{};

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
            {
                auto && [result, minimiser_count] = counter.count_query(seq, true);
                counts[counter_id] += result;

                size_t const threshold = arguments.threshold != 0.0 ?
                                            static_cast<size_t>(minimiser_count * arguments.threshold) :
                                            precomp_thresholds[std::min(minimiser_count < min_number_of_minimisers ?
                                                                            0 :
                                                                            minimiser_count - min_number_of_minimisers,
                                                                        max_number_of_minimisers -
                                                                            min_number_of_minimisers)] + 2;

                result_string.clear();
                result_string += id;
                result_string += '\t';

                size_t current_bin{0};
                for (auto && count : counts[counter_id++])
                {
                    if (count > threshold)
                    {
                        result_string += std::to_string(current_bin);
                        result_string += ',';
                    }
                    ++current_bin;
                }
                result_string += '\n';
                synced_out.write(result_string);
            }
        };

        do_parallel(output_task, records.size(), arguments.threads, compute_time);
    }

    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "IBF I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << ibf_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }
}

template <bool compressed>
void run_program_single(cmd_arguments const & arguments)
{
    ibf_builder<compressed> builder{&arguments};
    auto tbd = builder.ibf();

    std::ifstream is{arguments.ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&] ()
    {
        auto start = std::chrono::high_resolution_clock::now();
        iarchive(tbd);
        auto end = std::chrono::high_resolution_clock::now();
        ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    sync_out synced_out{arguments.out_file};

    size_t const kmers_per_window = arguments.window_size - arguments.kmer_size + 1;
    size_t const kmers_per_pattern = arguments.pattern_size - arguments.kmer_size + 1;
    size_t const min_number_of_minimisers = kmers_per_window == 1 ? kmers_per_pattern :
                                                std::ceil(kmers_per_pattern / static_cast<double>(kmers_per_window));
    size_t const max_number_of_minimisers = arguments.pattern_size - arguments.window_size + 1;
    std::vector<size_t> const precomp_thresholds = compute_simple_model(arguments);

    auto worker = [&] (size_t const start, size_t const end)
    {
        auto counter = tbd.template counting_agent<uint16_t>();
        std::string result_string{};

        for (auto && [id, seq] : records | seqan3::views::slice(start, end))
        {
            result_string.clear();
            result_string += id;
            result_string += '\t';

            auto && [result, minimiser_count] = counter.count_query(seq, true);
            size_t current_bin{0};

            size_t const threshold = arguments.threshold != 0.0 ?
                                         static_cast<size_t>(minimiser_count * arguments.threshold) :
                                         precomp_thresholds[std::min(minimiser_count < min_number_of_minimisers ?
                                                                         0 :
                                                                         minimiser_count - min_number_of_minimisers,
                                                                     max_number_of_minimisers -
                                                                         min_number_of_minimisers)] + 2;

            for (auto && count : result)
            {
                if (count > threshold)
                {
                    result_string += std::to_string(current_bin);
                    result_string += ',';
                }
                ++current_bin;
            }
            result_string += '\n';
            synced_out.write(result_string);
        }
    };

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        do_parallel(worker, records.size(), arguments.threads, compute_time);
    }

    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "IBF I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << ibf_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & arguments)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Search reads in a minimiser IBF.";
    parser.info.version = "1.0.0";
    parser.add_positional_option(arguments.query_file,
                                 "Please provide a path the FASTQ file.");
    parser.add_positional_option(arguments.ibf_file,
                                 "Please provide a valid path to an IBF. Parts: Without suffix _0");
    parser.add_option(arguments.out_file,
                      '\0',
                      "output",
                      "Please provide a valid path to the output.", seqan3::option_spec::DEFAULT);
    parser.add_option(arguments.window_size,
                      '\0',
                      "window",
                      "Choose the window size.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 35184372088832});
    parser.add_option(arguments.kmer_size,
                      '\0',
                      "kmer",
                      "Choose the kmer size.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(arguments.threads,
                      '\0',
                      "threads",
                      "Choose the number of threads.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 2048});
    parser.add_option(arguments.parts,
                      '\0',
                      "parts",
                      "Splits the IBF in this many parts. Must be a power of 2.");
    parser.add_option(arguments.errors,
                      '\0',
                      "error",
                      "Choose the number of errors.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{0, 5});
    parser.add_option(arguments.tau,
                      '\0',
                      "tau",
                      "Threshold for probabilistic models.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.threshold,
                      '\0',
                      "threshold",
                      "If set this threshold is used instead of the probabilistic models.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(arguments.pattern_size,
                      '\0',
                      "pattern",
                      "Choose the pattern size. Default: Use median of sequence lengths in query file.");
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed IBF.");
    parser.add_flag(arguments.write_time,
                    '\0',
                    "time",
                    "Write timing file.",
                    seqan3::option_spec::ADVANCED);
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"search_min_ibf", argc, argv, false};
    cmd_arguments arguments{};
    initialize_argument_parser(myparser, arguments);
    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    if (arguments.kmer_size > arguments.window_size)
        throw seqan3::argument_parser_error{"The kmer size cannot be bigger than the window size."};

    if (!arguments.pattern_size)
    {
        std::vector<uint64_t> sequence_lengths{};
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> query_in{arguments.query_file};
        for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
        {
            sequence_lengths.push_back(std::ranges::size(seq));
        }
        std::sort(sequence_lengths.begin(), sequence_lengths.end());
        arguments.pattern_size = sequence_lengths[sequence_lengths.size()/2];
    }

    if (arguments.parts == 1)
    {
        if (arguments.compressed)
            run_program_single<true>(arguments);
        else
            run_program_single<false>(arguments);
    }
    else
    {
        if (arguments.compressed)
            run_program_multiple<true>(arguments);
        else
            run_program_multiple<false>(arguments);
    }

    return 0;
}
