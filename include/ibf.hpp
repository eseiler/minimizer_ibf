// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_ibf/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides binning_directory.
 */

// See binning_directory PR https://github.com/seqan/seqan3/pull/920
// See DREAM Index PR https://github.com/seqan/seqan3/pull/851

#pragma once

#include <sdsl/bit_vectors.hpp>

#include <seqan/seq_io.h>

#include <strong_types.hpp>

/*!\brief The IBF binning directory.
 *
 * \details
 *
 * ### Binning Directory
 *
 * A binning directory is a data structure that can be used to determine set membership for elements.
 * For example, a common use case is dividing a database into a fixed number (e.g. 1024) bins by some means
 * of clustering (e.g. taxonomic binning or k-mer similarity clustering for genomic sequences).
 * For a query, the binning directory can now answer in which bins the query (probably) occurs.
 * In SeqAn we provide the Interleaved Bloom Filter (IBF) as underlying data structure.
 *
 * ### IBF
 *
 * The Interleaved Bloom Filter is a probabilistic data structure that extends the Bloom Filter.
 * A Bloom Filter can be thought of as a bitvector of length n and h hash functions and is used to determine set
 * membership. To insert data, the data is hashed by the h hash functions (returning values in [0, n)) and the
 * corresponding h positions in the bitvector are set to 1. To query data, i.e. to determine whether the query belongs
 * to the set the Bloom Filter was built for, the query is hashed by the same h hash functions and the corresponding
 * positions are checked. If all h positions contain a 1, the query is (probably) in the data set.
 * Since the Bloom Filter has variable length, the hashing is not bijective, i.e. it may return true for a set
 * membership query even though the query was never inserted into the Bloom Filter. Note that the Bloom Filter
 * will always return true if the query was inserted, i.e. there may be false positives, but no false negatives.
 *
 * The Interleaved Bloom Filter now applies the concept of a Bloom Filter to multiple sets and provides a *global*
 * data structure to determine set membership of a query in b data sets/bins.
 * Conceptually, a Bloom Filter is created for each bin using the same fixed length and fixed hash functions for each
 * filter. The resulting b Bloom Filters are then interleaved such that the i'th bit if each Bloom Filter are adjacent
 * to each other:
 * ```
 * Bloom Filter 0       Bloom Filter 1      Bloom Filter 2      Bloom Filter 3
 * |0.0|0.1|0.2|0.3|    |1.0|1.1|1.2|1.3|   |2.0|2.1|2.2|2.3|   |3.0|3.1|3.2|3.3|
 * ```
 * Where x.y denotes the y'th bit of the x'th Bloom Filter.
 * ```
 * Interleaved Bloom Filter
 * |0.0|1.0|2.0|3.0|0.1|1.1|2.1|3.1|0.2|1.2|2.2|3.2|0.3|1.3|2.3|3.3|
 * ```
 * A query can now be searched in all b bins by computing the h hash functions, retrieving the h sub-bitvectors of
 * length b starting at the positions indicated by the hash functions. The bitwise AND of these sub-bitvectors yields
 * the binningvector, a bitvector of length b where the i'th bit indicates set membership in the i'th bin.
 */
class binning_directory
{
private:
    // ===============================================================================================================
    // We fix some SeqAn2 related types since we won't use anything else for now.
    // ===============================================================================================================
    //!\brief The alphabet type.
    using alphabet_t = seqan::Dna;
    //!\brief The text type.
    using text_t = seqan::String<alphabet_t>;
    //!\brief The seqan::Shape type to compute rolling hashes with.
    using shape_t = seqan::Shape<alphabet_t, seqan::SimpleShape>;
    //!\brief Shape for computing the k-mers.
    shape_t hash_shape{};
    //!brief The k-mer size.
    size_t k{};

    //!\brief The number of bins.
    size_t nbins{};
    //!\brief The number of 64 bit integers used to represent the numerical value `nbins`.
    size_t bin_width{};
    //!\brief How big is nbins in multiples of 64.
    size_t block_size{};
    //!\brief How many blocks fit in the bitvector.
    size_t block_count{};
    //!\brief The number of hash functions.
    size_t num_hash{};
    //!\brief The bitvector.
    sdsl::bit_vector data{};
    //!\brief Precalculated values for hashing.
    std::vector<size_t> pre_hash{};
    //!\brief Shift value for hashing.
    static constexpr size_t shift{27};
    //!\brief Seed for hashing.
    static constexpr size_t seed{0x90b45d39fb6da1fa};

    /*!\brief Perturbs a value and fits it into the vector.
    * \param h The value to process.
    */
    inline constexpr void hash_and_fit(size_t & h) const
    {
        h ^= h >> shift; // Basically Fibonacci hashing (TODO optimize)
        h %= block_count;
        h *= block_size;
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    binning_directory() = default;                                      //!< Defaulted.
    binning_directory(binning_directory const &) = default;             //!< Defaulted.
    binning_directory & operator=(binning_directory const &) = default; //!< Defaulted.
    binning_directory(binning_directory &&) = default;                  //!< Defaulted.
    binning_directory & operator=(binning_directory &&) = default;      //!< Defaulted.
    ~binning_directory() = default;                                     //!< Defaulted.

    /*!\brief Construct using number of bins, bitvector size and the number of hash functions.
     * \param bins_     The number of bins.
     * \param bits_     The bitvector size.
     * \param num_hash_ The number of hash functions. Default 2.
     */
    binning_directory(bins bins_, bits bits_, kmer k_, hashes num_hash_ = hashes{2})
    {
        assert(bins_.v > 0);
        assert(bits_.v > 0);
        assert(num_hash_.v > 0);
        assert(k_.v > 0);
        if (bits_.v < bins_.v)
            throw std::logic_error{"There must be at least one bit per bin. Increase the current size of " +
                                   std::to_string(bits_.v) + " to at least " + std::to_string(bins_.v) + "."};
        nbins = bins_.v;
        k = k_.v;
        num_hash = num_hash_.v;
        data = sdsl::bit_vector(bits_.v);
        bin_width   = (nbins + 63) >> 6;
        block_size  = bin_width << 6;
        block_count = bits_.v / block_size;

        pre_hash.resize(num_hash);

        for (size_t i = 0; i < num_hash; ++i)
            pre_hash[i] = i ^ (k * seed);  // In the SeqAn2 implementation, we use `k`.
                                           // In SeqAn3 the IBF is decoupled from the hash value computation, so we
                                           // use a fixed value there (15).
    }
    //!\}

    /*!\brief Inserts a value into a specific bin.
     * \param h   The raw hash value to process.
     * \param bin The bin to insert into.
     */
    inline void set(size_t const h, size_t const bin)
    {
        assert(bin < nbins);
        // Sets respective positions in bitvector TODO SIMD?
        for (size_t const val : pre_hash)
        {
            size_t idx = val * h;
            hash_and_fit(idx);
            idx += bin;
            assert (idx < data.size());
            data[idx] = 1;
        }
    }

    /*!\brief Sets the number of bins.
     * \param new_bins The new number of bins.
     * \throws std::invalid_argument If passed number of bins is smaller than current number of bins.
     * \attention The new number of new bins must be greater or equal to the current number of bins.
     *
     * The resulting binning_directory has an increased size proportional to the increase in the bin_width, e.g.
     * resizing a binning_directory with 40 bins to 73 bins also increases the bin_width from 64 to 128 and hence the
     * new binning_directory will be twice the size.
     * This increase in size is necessary to avoid invalidating all computed hash functions.
     * If you want to add more bins while keeping the size constant, you need to rebuild the binning_directory.
     */
    inline void resize(size_t const new_bins)
    {
        if (new_bins < nbins)
        {
            throw std::invalid_argument{"The new number of new bins must be greater"
                                        " or equal to the current number of bins."};
        }

        size_t bin_width_ = (new_bins + 63) >> 6;

        if (bin_width_ == bin_width)
        {
            nbins = new_bins;
            return;
        }

        size_t block_size_ = bin_width_ << 6;
        size_t bits_       = block_count * block_size_;

        size_t idx_{bits_}, idx{data.size()};
        size_t delta = block_size_ - block_size + 64;

        data.resize(bits_);

        for (size_t i = idx_, j = idx; j > 0; i -= block_size_, j -= block_size)
        {
            size_t stop = i - block_size_;

            for (size_t ii = i - delta, jj = j - 64; stop && ii >= stop; ii -= 64, jj -= 64)
            {
                uint64_t old = data.get_int(jj);
                data.set_int(jj, 0);
                data.set_int(ii, old);
            }
        }

        nbins = new_bins;
        bin_width = bin_width_;
        block_size = block_size_;
    }

    /*!\brief Returns the number of bins in the binning directory.
     * \returns The number of bins.
     */
    inline size_t get_bins() const
    {
        return nbins;
    }

    /*!\brief Returns the number of hash functions used in the binning directory.
     * \returns The number of hash functions.
     */
    inline size_t get_num_hash() const
    {
        return num_hash;
    }

    /*!\brief Determines set membership of a given value.
     * \param[in,out] result The sdsl::bit_vector to store the result into.
     * \param[in]     h      The raw hash value to process.
     * \returns A sdsl::bit_vector of size nbins where each position indicates the bin membership of the hash value.
     */
    [[nodiscard]] sdsl::bit_vector get(size_t const h, sdsl::bit_vector && result) const
    {
        std::vector<size_t> idx = pre_hash;

        for (size_t & i : idx)
        {
            i *= h;
            hash_and_fit(i);
        }

        // SIMD?!! std::span?
        for (size_t batch = 0; batch < bin_width; ++batch)
        {
           assert(idx[0] < data.size());
           size_t tmp = data.get_int(idx[0]);
           idx[0] += 64;

           for (size_t i = 1; i < num_hash; ++i)
           {
               assert(idx[i] < data.size());
               tmp &= data.get_int(idx[i]);
               idx[i] += 64;
           }

           result.set_int(batch<<6, tmp);
        }

        return std::move(result);
    }

    //!\overload
    [[nodiscard]] sdsl::bit_vector get(size_t const h) const
    {
        sdsl::bit_vector result(nbins);

        return get(h, std::move(result));
    }

    /*!\brief Returns the size of the underlying bitvector.
     * \returns The size in bits of the underlying bitvector.
     */
    inline size_t size() const noexcept
    {
        return data.size();
    }

    /*!\brief Test for equality.
     * \param lhs A binning_directory.
     * \param rhs binning_directory to compare to.
     * \returns `true` if equal, `false` otherwise.
     */
    friend bool operator==(binning_directory const & lhs, binning_directory const & rhs) noexcept
    {
        return std::tie(lhs.nbins, lhs.bin_width, lhs.block_size, lhs.block_count, lhs.num_hash, lhs.data) ==
               std::tie(rhs.nbins, rhs.bin_width, rhs.block_size, rhs.block_count, rhs.num_hash, rhs.data);
    }

    /*!\brief Test for inequality.
     * \param lhs A binning_directory.
     * \param rhs binning_directory to compare to.
     * \returns `true` if unequal, `false` otherwise.
     */
    friend bool operator!=(binning_directory const & lhs, binning_directory const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    // ===============================================================================================================
    // Additional functions that won't be part of the seqan3::binning_directory, but part of a higher level interface.
    // ===============================================================================================================

    /*!\brief Insert a text into a specific bin.
     * \param bin  The bin to insert the data into.
     * \param text The text to process.
     */
    void insert_data(size_t const bin, text_t const & text) noexcept
    {
        seqan::resize(hash_shape, k);
        seqan::hashInit(hash_shape, seqan::begin(text));
        for (size_t i = 0; i < seqan::length(text) - k + 1; ++i)
            set(seqan::hashNext(hash_shape, seqan::begin(text) + i), bin);
    }

    /*!\brief Count the k-mers of a query in all bins.
     * \param query The query to count the k-mers for.
     * \returns A std::vector<size_t> of size bins where each element is the k-mer count for the respecitve bin.
     */
    std::vector<size_t> count(text_t const & query) noexcept
    {
        std::vector<size_t> result(nbins, 0);
        sdsl::bit_vector tmp_vec(nbins);

        seqan::resize(hash_shape, k);
        seqan::hashInit(hash_shape, seqan::begin(query));

        for (size_t i = 0; i < seqan::length(query) - k + 1; ++i)
        {
            tmp_vec = get(seqan::hashNext(hash_shape, seqan::begin(query) + i), std::move(tmp_vec));

            // TODO SIMD
            size_t bin{0};
            for (size_t batch = 0; batch < bin_width; ++batch)
            {
                size_t tmp = tmp_vec.get_int(batch * 64);
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
        return result;
    }

};
