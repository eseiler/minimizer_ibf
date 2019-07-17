#include <iostream>
#include <ibf.hpp>

int main()
{
    seqan::DnaString text0{"ACGTCGACGTTTAG"};
    seqan::DnaString text1{"TGCAAACGGCTTCA"};

    // Create a binning_directory with 64 bins and 8192 bits.
    binning_directory ibf{bins{64}, bits{8192}, kmer{4}};

    // Insert 4-mers of text0 into bin 0.
    ibf.insert_data(0, text0);

    // Insert 4-mers of text1 into bin 33.
    ibf.insert_data(33, text1);

    // Count 4-mers of text0 in the binning_directory.
    size_t i{0};
    for (auto && x : ibf.count(text0))
    {
        if (x > 0)
            std::cout << "Bin " << i << " : " << x << '\n'; // Bin 0 : 11
        ++i;
    }

    // Count 4-mers of text1 in the binning_directory.
    i = 0;
    for (auto && x : ibf.count(text1))
    {
        if (x > 0)
            std::cout << "Bin " << i << " : " << x << '\n'; // Bin 33 : 11
        ++i;
    }

}
