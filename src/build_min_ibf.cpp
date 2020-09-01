#include <seqan3/std/charconv>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct cmd_arguments
{
    std::vector<std::filesystem::path> bin_path{};
    std::filesystem::path out_path{"./"};
    std::string size{};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    uint32_t w{23};
    uint8_t k{20};
    uint8_t parts{1u};
    uint8_t threads{1};
    bool binary{false};
    bool bz2{false};
    bool compressed{false};
    bool gz{false};
    bool precomputed{false};
};

template <bool compressed>
struct ibf_builder
{
    cmd_arguments const * const arguments;

    template <typename view_t = std::ranges::empty_view<int>>
    auto construct(view_t && restrict_view = std::ranges::empty_view<int>())
    {
        std::string extension{".fasta"};
        if (arguments->gz)
            extension += ".gz";
        if (arguments->bz2)
            extension += ".bz2";

        using sequence_file_t = seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>>;

        std::vector<sequence_file_t> technical_bins;
        technical_bins.reserve(arguments->bins);

        for (size_t i = 0; i < arguments->bins; ++i)
            technical_bins.emplace_back(arguments->bin_path[0] / ("bin_" + std::to_string(i) + extension));

        seqan3::ibf_config cfg{seqan3::bin_count{arguments->bins},
                               seqan3::bin_size{arguments->bits / arguments->parts},
                               seqan3::hash_function_count{arguments->hash},
                               arguments->threads};

        if constexpr (std::same_as<view_t, std::ranges::empty_view<int>>)
        {
            return seqan3::technical_binning_directory{std::move(technical_bins),
                                                       seqan3::views::minimiser_hash(seqan3::ungapped{arguments->k},
                                                                                     seqan3::window_size{arguments->w},
                                                                                     seqan3::seed{adjust_seed(arguments->k)}),
                                                       cfg};
        }
        else
        {
            return seqan3::technical_binning_directory{std::move(technical_bins),
                                                       seqan3::views::minimiser_hash(seqan3::ungapped{arguments->k},
                                                                                     seqan3::window_size{arguments->w},
                                                                                     seqan3::seed{adjust_seed(arguments->k)})
                                                           | restrict_view,
                                                       cfg};
        }
    }

    template <typename view_t = std::ranges::empty_view<int>>
        requires !compressed
    auto ibf(view_t && restrict_view = std::ranges::empty_view<int>())
    {
        assert(arguments != nullptr);
        return construct(std::move(restrict_view));
    }

    template <typename view_t = std::ranges::empty_view<int>>
        requires compressed
    auto ibf(view_t && restrict_view = std::ranges::empty_view<int>())
    {
        assert(arguments != nullptr);
        auto tmp = construct(std::move(restrict_view));

        return seqan3::technical_binning_directory<seqan3::data_layout::compressed,
                                                   typename decltype(tmp)::hash_adaptor_t,
                                                   seqan3::dna4>{std::move(tmp)};
    }
};

template <bool compressed>
void run_program(cmd_arguments const & arguments)
{
    ibf_builder<compressed> generator{&arguments};

    if (arguments.parts == 1u)
    {
        auto tbd = generator.ibf();
        std::ofstream os{arguments.out_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(tbd);
    }
    else
    {
        std::vector<std::vector<size_t>> association(arguments.parts);
        size_t next_power_of_four{4u};

        if (arguments.parts == 4u) // one-to-one
        {
            for (size_t i : std::views::iota(0u, arguments.parts))
                association[i] = std::vector<size_t>{i};
        }
        else if (arguments.parts == 2u) // More than 1 prefix per part
        {
            association[0] = std::vector<size_t>{0, 1};
            association[1] = std::vector<size_t>{2, 3};
        }
        else // Multiple prefixes per part
        {
            // How long must the suffix be such that 4^suffix_length >= arguments.parts
            size_t suffix_length{0};
            for (; 0b100 << (2 * suffix_length) < arguments.parts; ++suffix_length) {}
            next_power_of_four = 0b100 << (2 * suffix_length);

            size_t const prefixes_per_part = next_power_of_four / arguments.parts;

            for (size_t i : std::views::iota(0u, next_power_of_four))
                association[i/prefixes_per_part].push_back(i);
        }

        for (size_t part : std::views::iota(0u, arguments.parts))
        {
            size_t const mask{next_power_of_four - 1};
            auto filter_view = std::views::filter([&] (auto && hash)
                { return std::ranges::find(association[part], hash & mask) != association[part].end(); });

            auto tbd = generator.ibf(filter_view);
            std::filesystem::path out_path{arguments.out_path};
            out_path += "_" + std::to_string(part);
            std::ofstream os{out_path, std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(tbd);
        }
    }
}

inline void compute_minimisers(cmd_arguments const & arguments)
{
    auto minimiser_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.k},
                                                        seqan3::window_size{arguments.w},
                                                        seqan3::seed{adjust_seed(arguments.k)});
    std::unordered_map<uint64_t, uint8_t> hash_table{}; // storage for minimisers
    uint64_t count{0};
    uint64_t filesize{0};
    uint16_t const default_cutoff{50};
    uint16_t cutoff{default_cutoff};
    // Cutoffs and bounds from Mantis
    // Mantis ignores k-mers which appear less than a certain cutoff. The cutoff is based on the file size of a
    // fastq gzipped file. So small files have only a cutoff of 1 while big files have a cutoff value of 50.
    std::vector<uint16_t> const cutoffs{1,3,10,20};
    std::vector<uint64_t> const cutoff_bounds{314572800, 524288000, 1073741824, 3221225472};

    for (auto const & file : arguments.bin_path)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{file};

        for (auto & [seq] : fin)
            for (auto && hash : seq | minimiser_view)
                hash_table[hash] = std::min<uint8_t>(254u, hash_table[hash] + 1);
                // The hash table stores how often a minimiser appears. It does not matter whether a minimiser appears
                // 50 times or 2000 times, it is stored regardless because the biggest cutoff value is 50. Hence
                // the hash table stores only values up until 254 to save memory.

        // The filesize is multiplied by two, because Mantis' filesize is based on fastq files, while we use fasta files
        // which are smaller because the quality information is missing. A multiplication of 2 should return roughly
        // the size of the fastq file. Note: Because Mantis cutoffs are based on gzipped files, this estimation makes
        // only sense when gzipped fasta files are used.
        filesize = std::filesystem::file_size(file) * 2;

        for (size_t k = 0; k < cutoff_bounds.size(); ++k)
        {
            if (filesize <= cutoff_bounds[k])
            {
                cutoff = cutoffs[k];
                break;
            }
        }

        // Store binary file
        std::ofstream outfile{arguments.out_path.string() + file.stem().string() + ".minimiser", std::ios::binary};
        for (auto && hash : hash_table)
        {
            if (hash.second > cutoff)
            {
                outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
                ++count;
            }
        }

        // Store header file
        std::ofstream headerfile{arguments.out_path.string() + file.stem().string() + ".header"};
        headerfile << static_cast<uint64_t>(arguments.k) << '\t'
                   << arguments.w << '\t'
                   << cutoff << '\t'
                   << count << '\n';

        count = 0;
        cutoff = default_cutoff;
        hash_table.clear();
    }
}

template <bool compressed>
void build_from_binary(cmd_arguments const & arguments)
{
    seqan3::ibf_config cfg{seqan3::bin_count{arguments.bins},
                           seqan3::bin_size{arguments.bits / arguments.parts},
                           seqan3::hash_function_count{arguments.hash},
                           arguments.threads};

    seqan3::technical_binning_directory tbd{std::vector<std::vector<seqan3::dna4>>{},
                                            seqan3::views::minimiser_hash(seqan3::ungapped{arguments.k},
                                                                          seqan3::window_size{arguments.w},
                                                                          seqan3::seed{adjust_seed(arguments.k)}),
                                            cfg,
                                            true};

    uint64_t num;

    for (uint64_t cur_bin = 0; cur_bin < arguments.bin_path.size(); ++cur_bin)
    {
        std::ifstream infile{arguments.bin_path[cur_bin], std::ios::binary};

        while(infile.read(reinterpret_cast<char*>(&num), sizeof(num)))
            tbd.emplace(num, seqan3::bin_index{cur_bin});
    }

    if constexpr (compressed)
    {
        seqan3::technical_binning_directory<seqan3::data_layout::compressed,
                                            typename decltype(tbd)::hash_adaptor_t,
                                            seqan3::dna4> ctbd{std::move(tbd)};

        std::ofstream os{arguments.out_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(ctbd);
    }
    else
    {
        std::ofstream os{arguments.out_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(tbd);
    }
}

void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & arguments)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Build an Interleaved Bloom Filter using minimisers.";
    parser.info.version = "0.0.1";
    parser.add_option(arguments.bin_path,
                      '\0',
                      "input",
                      "Please provide a path to a directory containing one FASTA file for each bin. Or provide a list "
                          "of binary minimiser files as output by build_min_ibf --binary.",
                      seqan3::option_spec::REQUIRED);
    parser.add_option(arguments.out_path,
                      '\0',
                      "output",
                      "Please provide a valid output path. Acts as a prefix for -binary. Defaults to \"ibf.out\" when "
                          "not using --binary.",
                      seqan3::option_spec::DEFAULT);
    parser.add_option(arguments.w,
                      '\0',
                      "window",
                      "Choose the window size.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 35184372088832});
    parser.add_option(arguments.k,
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
    parser.add_option(arguments.bins,
                      '\0',
                      "bins",
                      "Choose the number of bins.",
                      seqan3::option_spec::REQUIRED,
                      seqan3::arithmetic_range_validator{1, 65536});
    parser.add_option(arguments.bits,
                      '\0',
                      "bits",
                      "Choose the size in bits of one bin. Mutually exclusive with --size.",
                      seqan3::option_spec::DEFAULT,
                                      seqan3::arithmetic_range_validator{1, 35184372088832});
    parser.add_option(arguments.size,
                      '\0',
                      "size",
                      "Choose the size of the resulting IBF. Mutually exclusive with --bits. "
                          "Allowed suffices: {k, m, g, t}, e.g., --size 8g.");
    parser.add_option(arguments.parts,
                      '\0',
                      "parts",
                      "Splits the IBF in this many parts. Must be a power of 2.");
    parser.add_option(arguments.hash,
                      '\0',
                      "hash",
                      "Choose the number of hashes.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 4});
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed IBF.");
    parser.add_flag(arguments.gz,
                    '\0',
                    "gz",
                    "Expect FASTA files to be gz compressed.");
    parser.add_flag(arguments.bz2,
                    '\0',
                    "bz2",
                    "Expect FASTA files to be bz2 compressed.");
    parser.add_flag(arguments.binary,
                    '\0',
                    "binary",
                    "Only precompute minimisers. Does not create the IBF.");
    parser.add_flag(arguments.precomputed,
                    '\0',
                    "precomputed",
                    "Use precomputed minimisers to create IBF.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"build_ibf", argc, argv, false};
    cmd_arguments arguments{};
    initialise_argument_parser(myparser, arguments);
    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    // ==========================================
    // Set default for output path.
    // ==========================================
    if (!arguments.binary && arguments.out_path == "./")
        arguments.out_path = "./ibf.out";

    // ==========================================
    // Various checks.
    // ==========================================
    if (arguments.gz && arguments.bz2)
        throw seqan3::argument_parser_error{"Files cannot be both gz and bz2 compressed."};

    if (!arguments.precomputed && arguments.bin_path.size() > 1)
    {
        throw seqan3::argument_parser_error{"Too many paths for --input. Only provide one path to a directory when not "
                                            "using --precomputed."};
    }

    if (arguments.precomputed && arguments.bin_path.size() > arguments.bins)
        throw seqan3::argument_parser_error{"There are more paths given via --input than there are bins."};

    if (arguments.k > arguments.w)
        throw seqan3::argument_parser_error{"The kmer size cannot be bigger than the window size."};

    // ==========================================
    // Validate input files.
    // ==========================================
    {
        if (arguments.precomputed)
        {
            for (auto const & file_path : arguments.bin_path)
            {
                std::error_code ec{};
                std::filesystem::exists(file_path, ec);

                if (ec)
                {
                    throw seqan3::argument_parser_error{"The file " + file_path.string() +
                                                        "could not be opened/read.\n" + ec.message()};
                }
            }
        }
        else
        {
            std::string extension{};

            if (arguments.gz)
                extension += ".gz";

            if (arguments.bz2)
                extension += ".bz2";

            for (size_t i = 0; i < arguments.bins; ++i)
            {
                std::filesystem::path file_path{arguments.bin_path[0] / ("bin_" + std::to_string(i) + extension)};
                std::error_code ec{};
                std::filesystem::exists(file_path, ec);

                if (ec)
                {
                    throw seqan3::argument_parser_error{"The file " + file_path.string() +
                                                        "could not be opened/read.\n" + ec.message()};
                }
            }
        }
    }

    // ==========================================
    // Process --size.
    // ==========================================
    arguments.size.erase(std::remove(arguments.size.begin(), arguments.size.end(), ' '), arguments.size.end());

    // Probably not default. https://github.com/seqan/seqan3/pull/1859
    if (arguments.bits != 4096u && !arguments.size.empty())
        throw seqan3::argument_parser_error{"Either set --bits or --size."};

    if (!arguments.size.empty())
    {
        size_t multiplier{};

        switch (std::tolower(arguments.size.back()))
        {
            case 't':
                multiplier = 8ull * 1024ull * 1024ull * 1024ull * 1024ull;
                break;
            case 'g':
                multiplier = 8ull * 1024ull * 1024ull * 1024ull;
                break;
            case 'm':
                multiplier = 8ull * 1024ull * 1024ull;
                break;
            case 'k':
                multiplier = 8ull * 1024ull;
                break;
            default:
                throw seqan3::argument_parser_error{"Use {k, m, g, t} to pass size. E.g., --size 8g."};
        }

        size_t size{};
        std::from_chars(arguments.size.data(), arguments.size.data() + arguments.size.size() - 1, size);
        size *= multiplier;
        arguments.bits = size / (((arguments.bins + 63) >> 6) << 6);
    }

    // ==========================================
    // Run correct method.
    // ==========================================
    if (arguments.binary)
    {
        compute_minimisers(arguments);
        return 0;
    }
    if (arguments.precomputed)
    {
        if (arguments.compressed)
            build_from_binary<true>(arguments);
        else
            build_from_binary<false>(arguments);
        return 0;
    }
    if (arguments.compressed)
        run_program<true>(arguments);
    else
        run_program<false>(arguments);
    return 0;
}
