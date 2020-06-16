#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <iostream>
#include <minimizer.hpp>

struct cmd_arguments
{
    std::vector<std::filesystem::path> bin_path{};
    std::filesystem::path out_path{"./"};
    uint64_t w{23};
    uint8_t k{20};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    bool gz{false};
    bool bz2{false};
    bool binary{false};
    bool precalculated{false};
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

void fill_ibf(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> & ibf,
              std::vector<std::filesystem::path> const & dir_path,
              uint8_t const n_k,
              uint64_t const n_w,
              std::string const & extension)
{
    minimizer mini{window{n_w}, kmer{n_k}};

    for (uint64_t cur_bin = 0; cur_bin < ibf.bin_count(); ++cur_bin)
    {
        std::filesystem::path bin_path{dir_path[0]};
        bin_path /= ("bin_" + std::to_string(cur_bin) + extension);

        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{bin_path};

        for (auto & [seq] : fin)
        {
            mini.compute(seq);
            for (auto && hash : mini.minimizer_hash)
                ibf.emplace(hash, seqan3::bin_index{cur_bin});
        }
    }
}

// Reads precomputed minimisers from binary files.
void fill_ibf(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> & ibf,
              std::vector<std::filesystem::path> const & in_path)
{
    std::ifstream infile;

    for (uint64_t cur_bin = 0; cur_bin < in_path.size(); ++cur_bin)
    {
        infile.open(in_path[cur_bin], std::ios::binary);

        if (!infile.is_open())
        {
            throw std::filesystem::filesystem_error("File not found!",
                                                    in_path[cur_bin],
                                                    std::make_error_code(std::errc::bad_file_descriptor));
        }

        uint64_t num;
        while(infile.read(reinterpret_cast<char*>(&num), sizeof(num)))
            ibf.emplace(num, seqan3::bin_index{cur_bin});

        infile.close();
    }
}

void run_program(cmd_arguments & args)
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{args.bins},
                                         seqan3::bin_size{args.bits},
                                         seqan3::hash_function_count{args.hash}};

    std::string extension{".fasta"};
    if (args.gz)
        extension += ".gz";
    if (args.bz2)
        extension += ".bz2";

    if (args.precalculated)
        fill_ibf(ibf, args.bin_path);
    else
        fill_ibf(ibf, args.bin_path, args.k, args.w, extension);

    if (args.out_path.string() == "./")
        args.out_path = std::filesystem::path{"./out.ibf"};

    std::ofstream os{args.out_path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}

void get_minimizer_files(std::vector<std::filesystem::path> const & in_paths,
                         std::filesystem::path const & out_path,
                         uint8_t const n_k,
                         uint64_t const n_w)
{
    std::unordered_map<uint64_t, uint8_t> hash_table{}; // storage for minimizers
    std::ofstream outfile;
    std::ofstream headerfile;
    uint64_t count{0};
    uint64_t filesize{0};
    uint16_t cutoff{50}; // Default cutoff.
    // Cutoffs and bounds from Mantis
    // Mantis ignores k-mers which appear less than a certain cutoff. The cutoff is based on the file size of a
    // fastq gzipped file. So small files have only a cutoff of 1, while big files have a cutoff value of 50.
    std::vector<uint16_t> cutoffs{1,3,10,20};
    std::vector<uint64_t> cutoff_bounds{314572800, 524288000, 1073741824, 3221225472};

    minimizer mini{window{n_w}, kmer{n_k}};

    for (auto & file : in_paths)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{file};

        for (auto & [seq] : fin)
        {
            mini.compute(seq);
            for (auto && hash : mini.minimizer_hash)
                // The hash table stores how often a minimizer appears. It does not matter if a minimizer appears
                // 50 times or 2000 times, it is stored regardless because the biggest cutoff value is 50, therefore
                // the hash table stores only values up until 254 to save memory.
                hash_table[hash] = std::min<uint8_t>(254u, hash_table[hash] + 1);
        }

        // The filesize is multiplied by two, because Mantis' filesize is based on fastq files, while we use fasta files,
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
        outfile.open(std::string(out_path) + std::string(file.stem()) + ".minimizer", std::ios::binary);
        for (auto & hash : hash_table)
        {
            if (hash.second > cutoff)
            {
                outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
                ++count;
            }
        }
        outfile.close();

        // Store header file
        headerfile.open(std::string(out_path) + std::string(file.stem()) + ".header");
        headerfile << static_cast<uint64_t>(n_k) << "\t" << n_w << "\t" << cutoff << "\t" << count << "\n";
        headerfile.close();

        count = 0;
        cutoff = 50;
        hash_table.clear();
    }
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Build an IBF using minimizers.";
    parser.info.version = "1.0.0";
    parser.add_positional_option(args.bin_path, "Please provide a path to a directory containing one FASTA file for "
                                                "each bin. Or provide a list of binary minimizer files as output by "
                                                "build_min_ibf --binary.");
    parser.add_option(args.out_path, 'o', "out", "Please provide a valid output path.", seqan3::option_spec::DEFAULT);
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
    parser.add_flag(args.binary, '\0', "binary", "Only precompute minimizers. Does not create the IBF.");
    parser.add_flag(args.precalculated, '\0', "pre", "Use precomputed minimizers to create IBF.");
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

    if (args.binary)
    {
        get_minimizer_files(args.bin_path, args.out_path, args.k, args.w);
        return 0;
    }

    if (args.gz && args.bz2)
        throw seqan3::argument_parser_error{"Files cannot be both gz and bz2 compressed."};

    if (args.k > args.w)
        throw seqan3::argument_parser_error{"The kmer size cannot be bigger than the window size."};

    run_program(args);
    return 0;
}
