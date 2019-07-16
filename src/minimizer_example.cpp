#include <iostream>
#include <minimizer.hpp>

int main()
{
    seqan::DnaString test{"ACGTCGACGTTTAG"};

    {
    // With xor. (Default)
    minimizer h{window{8}, kmer{4}}; // Same as minimizer<use_xor::yes> h{window{8}, kmer{4}};
    h.compute(test);
    for (auto && x : h.minimizer_begin)
        std::cout << x << '\t'; // 3 5 9
    std::cout << '\n';
    for (auto && x : h.minimizer_end)
        std::cout << x << '\t'; // 6 8 12
    std::cout << '\n';
    for (auto && x : h.minimizer_hash)
        std::cout << x << '\t'; // 10322096095657499142 10322096095657499224 10322096095657499166
    std::cout << '\n';
    }
    std::cout << '\n';
    {
    // Without xor.
    minimizer<use_xor::no> h{window{8}, kmer{4}};
    h.compute(test);
    for (auto && x : h.minimizer_begin)
        std::cout << x << '\t'; // 0 2 6 7 8
    std::cout << '\n';
    for (auto && x : h.minimizer_end)
        std::cout << x << '\t'; // 3 5 9 10 11
    std::cout << '\n';
    for (auto && x : h.minimizer_hash)
        std::cout << x << '\t'; // 27 97 27 6 1
    std::cout << '\n';
    }
    std::cout << '\n';
    {
    // Custom seed.
    // With our without `<use_xor::yes>`
    minimizer h{window{8}, kmer{4}, 0xFFFF'FFFF'FFFF'FFFFULL};
    h.compute(test);
    for (auto && x : h.minimizer_begin)
        std::cout << x << '\t'; // 3 8 9
    std::cout << '\n';
    for (auto && x : h.minimizer_end)
        std::cout << x << '\t'; // 6 11 12
    std::cout << '\n';
    for (auto && x : h.minimizer_hash)
        std::cout << x << '\t'; // 18446744073709551399	18446744073709551424	18446744073709551363
    std::cout << '\n';
    }
}
