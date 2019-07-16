#include <benchmark/benchmark.h>

#include <minimizer.hpp>

seqan::CharString example_file{std::string{{BASE_DIR}} + "/data/example.fasta"};

template <use_xor xor_t>
static void compute_minimizer(benchmark::State & state)
{
    uint64_t w = static_cast<uint64_t>(state.range(0));
    uint64_t k = static_cast<uint64_t>(state.range(1));
    double io_time{0.0};
    double hash_time{0.0};

    minimizer<xor_t> mini{window{w}, kmer{k}};

    for (auto _ : state)
    {
        seqan::CharString id;
        seqan::String<seqan::Dna> seq;
        seqan::SeqFileIn seqFileIn;
        if (!seqan::open(seqFileIn, seqan::toCString(example_file)))
        {
            seqan::CharString msg = "Unable to open contigs file: ";
            seqan::append(msg, seqan::CharString(example_file));
            throw seqan::toCString(msg);
        }
        while(!seqan::atEnd(seqFileIn))
        {
            auto start = std::chrono::high_resolution_clock::now();
            seqan::readRecord(id, seq, seqFileIn);
            auto end = std::chrono::high_resolution_clock::now();
            io_time += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            mini.compute(seq);
            end = std::chrono::high_resolution_clock::now();
            hash_time += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
        }
    }
    state.counters["io_time"] = io_time;
    state.counters["hash_time"] = hash_time;
}

static void minimizer_arguments(benchmark::internal::Benchmark* b)
{
    for (int32_t w = 20; w < 23; ++w)
    {
        for (int32_t k = w - 3; k <= w; ++k)
        {
            b->Args({w, k});
        }
    }
}

BENCHMARK_TEMPLATE(compute_minimizer, use_xor::no)->Apply(minimizer_arguments);
BENCHMARK_TEMPLATE(compute_minimizer, use_xor::yes)->Apply(minimizer_arguments);

BENCHMARK_MAIN();
