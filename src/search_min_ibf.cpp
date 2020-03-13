#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/async_input_buffer.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <iostream>
#include <minimizer.hpp>


struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

void run_program(std::filesystem::path const & query_path,
                 std::filesystem::path const & ibf_path,
                 uint8_t const n_k,
                 uint64_t const n_w,
                 uint8_t const)
{
    seqan3::interleaved_bloom_filter ibf;

    std::ifstream is{ibf_path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{query_path};

    // create the async buffer around the input file
    // spawns a background thread that tries to keep eight records in the buffer
    auto v = fin | seqan3::views::async_input_buffer(8);

    // create a lambda function that iterates over the async buffer when called
    // (the buffer gets dynamically refilled as soon as possible)
    minimizer<use_xor::yes> mini{window{n_w}, kmer{n_k}};
    std::vector<size_t> result(ibf.bin_count(), 0);

    for (auto & [seq] : v)
    {
        result.clear();
        mini.compute(seq);

        for (auto && hash : mini.minimizer_hash)
        {
            auto & res = ibf.bulk_contains(hash);

            size_t bin{0};
            for (size_t batch = 0; batch < ((ibf.bin_count() + 63) >> 6); ++batch)
            {
                size_t tmp = res.get_int(batch * 64);
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
    }
}

struct cmd_arguments
{
    std::filesystem::path query_path{};
    std::filesystem::path ibf_path{};
    uint64_t w{23};
    uint8_t k{20};
    uint8_t t{1};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Search reads in a minimizer IBF.";
    parser.info.version = "1.0.0";
    parser.add_positional_option(args.query_path, "Please provide a path the FASTQ file.");
    parser.add_positional_option(args.ibf_path, "Please provide a valid path to a minimizer IBF.");
    parser.add_option(args.w, '\0', "window", "Choose the window size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 1000});
    parser.add_option(args.k, '\0', "kmer", "Choose the kmer size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(args.t, '\0', "threads", "Choose the number of threads.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 2048});
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

    run_program(args.query_path, args.ibf_path, args.k, args.w, args.t);
    return 0;
}
