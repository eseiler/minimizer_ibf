#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <iostream>
#include <minimizer.hpp>


struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

void run_program(std::filesystem::path const & dir_path,
                 std::filesystem::path const & out_path,
                 uint8_t const n_k,
                 uint64_t const n_w,
                 uint64_t const n_bins,
                 uint64_t const n_bits,
                 uint64_t const n_hash)
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{n_bins},
                                         seqan3::bin_size{n_bits},
                                         seqan3::hash_function_count{n_hash}};
    minimizer mini{window{n_w}, kmer{n_k}};

    std::string extension{".fasta"};
    if (args.gz)
        extension += ".gz";
    if (args.bz2)
        extension += ".bz2";

    for (uint64_t cur_bin = 0; cur_bin < n_bins; ++cur_bin)
    {
        std::filesystem::path bin_path{dir_path};
        bin_path /= ("bin_" + std::to_string(cur_bin) + extension);

        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{bin_path};

        for (auto & [seq] : fin)
        {
            mini.compute(seq);
            for (auto && hash : mini.minimizer_hash)
                ibf.emplace(hash, seqan3::bin_index{cur_bin});
        }
    }

    std::ofstream os{out_path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}

struct cmd_arguments
{
    std::filesystem::path bin_path{};
    std::filesystem::path out_path{};
    uint64_t w{23};
    uint8_t k{20};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    bool gz{false};
    bool bz2{false};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Build an IBF using minimizers.";
    parser.info.version = "1.0.0";
    parser.add_positional_option(args.bin_path, "Please provide a path to a directory containing one FASTA file for each bin.");
    parser.add_positional_option(args.out_path, "Please provide a valid output path.");
    parser.add_option(args.w, '\0', "window", "Choose the window size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 1000});
    parser.add_option(args.k, '\0', "kmer", "Choose the kmer size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(args.bins, '\0', "bins", "Choose the number of bins.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 65536});
    parser.add_option(args.bits, '\0', "bits", "Choose the size in bits of one bin.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 35184372088832});
    parser.add_option(args.hash, '\0', "hash", "Choose the number of hashes.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 4});
    parser.add_flag(args.gz, '\0', "gz", "Expect FASTA files to be gz compressed.");
    parser.add_flag(args.bz2, '\0', "bz2", "Expect FASTA files to be bz2 compressed.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"build_min_ibf", argc, argv, false};
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

    if (args.gz && args.bz2)
        throw seqan3::argument_parser_error{"Files cannot be both gz and bz2 compressed."};

    run_program(args.bin_path, args.out_path, args.k, args.w, args.bins, args.bits, args.hash);
    return 0;
}
