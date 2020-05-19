#include <minimizer_model.hpp>

#define DEBUG 0

std::vector<size_t> compute_simple_model(cmd_arguments const & args)
{
    std::vector<size_t> precomp_thresholds;

    if (!do_cerealisation_in(precomp_thresholds, args))
    {
        precomp_thresholds = precompute_threshold(args.pattern_size,
                                                  args.window_size,
                                                  args.kmer_size,
                                                  args.errors,
                                                  args.tau);

        do_cerealisation_out(precomp_thresholds, args);
    }

    return precomp_thresholds;
}

void run_program(cmd_arguments const & args)
{
    seqan3::interleaved_bloom_filter ibf;

    std::ifstream is{args.ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{args.query_file};

#if DEBUG
    // Extract bin_number from `bin_[x]+`-format query file name, otherwise set to 0.
    size_t const bin_no = [&] () {
        try
        {
            return static_cast<size_t>(std::stoi(args.query_file.stem().string().substr(4)));
        }
        catch (std::exception const & e)
        {
            return static_cast<size_t>(0);
        }
    }();
#endif

    // create the async buffer around the input file
    // spawns a background thread that tries to keep eight records in the buffer
    auto sequence_input_buffer = fin | seqan3::views::async_input_buffer(8);

    // Some computations related to the minimizer threshold
    size_t const kmers_per_window = args.window_size - args.kmer_size + 1;
    size_t const kmers_per_pattern = args.pattern_size - args.kmer_size + 1;
    size_t const minimal_number_of_minimizers = kmers_per_window == 1 ? kmers_per_pattern :
                                                                        kmers_per_pattern / (kmers_per_window - 1);
    size_t const maximal_number_of_minimizers = args.pattern_size - args.window_size + 1;

    std::vector<size_t> const precomp_thresholds = compute_simple_model(args);

#if DEBUG
    // Counts the number of read sequences/queries
    std::atomic<size_t> seq_count{0};
    // Counts the number of hits in the searched bin; inferred from query file name
    std::atomic<size_t> hit_count{0};
#endif

    // create a lambda function that iterates over the async buffer when called
    // (the buffer gets dynamically refilled as soon as possible)
    auto worker = [&] ()
    {
        minimizer mini{window{args.window_size}, kmer{args.kmer_size}};
        decltype(ibf)::binning_bitvector result_buffer(ibf.bin_count());
        std::vector<size_t> result(ibf.bin_count(), 0);

        for (auto & [seq] : sequence_input_buffer)
        {
#if DEBUG
            ++seq_count;
#endif
            std::fill(result.begin(), result.end(), 0);
            mini.compute(seq);

            for (auto && hash : mini.minimizer_hash)
            {
                ibf.bulk_contains(hash, result_buffer);

                size_t bin{0};
                for (size_t batch = 0; batch < ((ibf.bin_count() + 63) >> 6); ++batch)
                {
                    size_t tmp = result_buffer.get_int(batch * 64);
                    if (tmp ^ (1ULL<<63))
                    {
                        while (tmp > 0)
                        {
                            uint8_t step = sdsl::bits::lo(tmp);
                            bin += step++;
                            tmp >>= step;
                            ++result[bin++];
                        }
                    }
                    else
                    {
                        ++result[bin + 63];
                    }
                }
            }

            size_t const minimizer_count = mini.minimizer_hash.size();
            size_t index = std::min(minimizer_count < minimal_number_of_minimizers ? 0 :
                                        minimizer_count - minimal_number_of_minimizers,
                                    maximal_number_of_minimizers - minimal_number_of_minimizers);
            auto threshold = precomp_thresholds[index];

            size_t current_bin{0};
            for (auto & count : result)
            {
                count = count >= threshold;
#if DEBUG
                    if (count && bin_no == current_bin)
                        ++hit_count;
#endif
                ++current_bin;
            }
        }
    };

    std::vector<decltype(std::async(std::launch::async, worker))> tasks;

    for (size_t i = 0; i < args.threads; ++i)
        tasks.emplace_back(std::async(std::launch::async, worker));

    for (auto && task : tasks)
        task.wait();

#if DEBUG
    std::cout << "seq_count " << seq_count.load() << '\t' << "hit_count " << hit_count.load() << '\n';
#endif
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Search reads in a minimizer IBF.";
    parser.info.version = "1.0.0";
    parser.add_positional_option(args.query_file, "Please provide a path the FASTQ file.");
    parser.add_positional_option(args.ibf_file, "Please provide a valid path to a minimizer IBF.");
    parser.add_option(args.window_size, '\0', "window", "Choose the window size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 1000});
    parser.add_option(args.kmer_size, '\0', "kmer", "Choose the kmer size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(args.threads, '\0', "threads", "Choose the number of threads.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 2048});
    parser.add_option(args.errors, '\0', "error", "Choose the number of errors.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{0, 5});
    parser.add_option(args.tau, '\0', "tau", "Threshold for probabilistic models.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{0, 1});
    parser.add_option(args.pattern_size, '\0', "pattern",
                      "Choose the pattern size. Default: Use median of sequence lengths in query file.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"search_min_ibf", argc, argv, false};
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);
    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    if (!args.pattern_size)
    {
        std::vector<uint64_t> sequence_lengths{};
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> query_in{args.query_file};
        for (auto & [seq] : query_in | seqan3::views::async_input_buffer(16))
        {
            sequence_lengths.push_back(std::ranges::size(seq));
        }
        std::sort(sequence_lengths.begin(), sequence_lengths.end());
        args.pattern_size = sequence_lengths[sequence_lengths.size()/2];
    }

    run_program(args);
    return 0;
}
