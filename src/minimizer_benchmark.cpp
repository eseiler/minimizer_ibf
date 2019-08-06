#include <benchmark/benchmark.h>

#include <minimizer.hpp>

seqan::CharString file{};

template <use_xor xor_t>
static void compute_minimizer(benchmark::State & state)
{
    uint64_t w = static_cast<uint64_t>(state.range(0));
    uint64_t k = static_cast<uint64_t>(state.range(1));
    std::chrono::duration<double, std::milli> io_time{0.0};
    std::chrono::duration<double, std::milli> hash_time{0.0};

    minimizer<xor_t> mini{window{w}, kmer{k}};

    for (auto _ : state)
    {
        seqan::CharString id;
        seqan::String<seqan::Dna> seq;
        seqan::SeqFileIn seqFileIn;
        if (!seqan::open(seqFileIn, seqan::toCString(file)))
        {
            seqan::CharString msg = "Unable to open contigs file: ";
            seqan::append(msg, seqan::CharString(file));
            throw std::runtime_error{seqan::toCString(msg)};
        }
        while(!seqan::atEnd(seqFileIn))
        {
            auto start = std::chrono::high_resolution_clock::now();
            seqan::readRecord(id, seq, seqFileIn);
            auto end = std::chrono::high_resolution_clock::now();
            io_time += end - start;

            start = std::chrono::high_resolution_clock::now();
            mini.compute(seq);
            end = std::chrono::high_resolution_clock::now();
            hash_time += end - start;
        }
    }
    state.counters["io_time_sum[ms]"] = io_time.count();
    state.counters["io_time_avg[ms]"] = io_time.count() / state.iterations();
    state.counters["hash_time_sum[ms]"] = hash_time.count();
    state.counters["hash_time_avg[ms]"] = hash_time.count() / state.iterations();
}

static void minimizer_arguments(benchmark::internal::Benchmark* b)
{
    for (int32_t w = 23; w <= 23; ++w)
    {
        for (int32_t k = w - 3; k <= w; ++k)
        {
            b->Args({w, k});
        }
    }
}

BENCHMARK_TEMPLATE(compute_minimizer, use_xor::no)->Apply(minimizer_arguments);
BENCHMARK_TEMPLATE(compute_minimizer, use_xor::yes)->Apply(minimizer_arguments);

int main(int argc, char** argv)
{
    benchmark::Initialize(&argc, argv);

    if (argc != 2)
        throw std::invalid_argument{"Please specify the input file."};

    file = seqan::CharString{argv[1]};
    --argc;

    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
