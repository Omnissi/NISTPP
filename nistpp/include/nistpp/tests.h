/**
 * @file
 *
 * @author Negodyaev Sergey (negodyaev.sergey@outlook.com)
 * @details This file contains declarations of all tests.
 */
#ifndef TESTS_H
#define TESTS_H

#include <cstddef>
#include "bits_storage.h"
#include "types.h"

namespace nistpp
{

/// @brief Frequency (Monobit) Test
/// @details The focus of the test is the proportion of zeroes and ones for the entire sequence. The purpose of this test
/// is to determine whether the number of ones and zeros in a sequence are approximately the same as would
/// be expected for a truly random sequence. The test assesses the closeness of the fraction of ones to 1/2, that
/// is, the number of ones and zeroes in a sequence should be about the same. All subsequent tests depend on
/// the passing of this test.
///
/// @param[in] data Class contained sequence for test.
/// @return @ref nistpp::return_t
return_t FrequencyTest(const BitsStorage& data);

/// @brief Frequency Test within a Block.
/// @details The focus of the test is the proportion of ones within M-bit blocks. The purpose of this test is to determine
/// whether the frequency of ones in an M-bit block is approximately M/2, as would be expected under an
/// assumption of randomness. For block size M=1, this test degenerates to test 1, the Frequency (Monobit) test.
///
/// @param[in] data Class contained sequence for test.
/// @param[in] M    The length of each block.
/// @return @ref nistpp::return_t
/// @throw std::invalid_argument If M > number of bits.
return_t BlockFrequencyTest(const BitsStorage& data, std::size_t M);

/// @brief Runs Test.
/// @details The focus of this test is the total number of runs in the sequence, where a run is an uninterrupted sequence
/// of identical bits. A run of length k consists of exactly k identical bits and is bounded before and after with
/// a bit of the opposite value. The purpose of the runs test is to determine whether the number of runs of
/// ones and zeros of various lengths is as expected for a random sequence. In particular, this test determines
/// whether the oscillation between such zeros and ones is too fast or too slow.
///
/// @param[in] data Class contained sequence for test.
/// @return @ref nistpp::return_t
return_t RunsTest(const BitsStorage& data);

/// @brief Test for the Longest Run of Ones in a Block.
/// @details The focus of the test is the longest run of ones within M-bit blocks. The purpose of this test is to
/// determine whether the length of the longest run of ones within the tested sequence is consistent with the
/// length of the longest run of ones that would be expected in a random sequence. Note that an irregularity in
/// the expected length of the longest run of ones implies that there is also an irregularity in the expected
/// length of the longest run of zeroes. Therefore, only a test for ones is necessary.
///
/// @param[in] data Class contained sequence for test.
/// @return @ref nistpp::return_t
return_t LongestRunOfOnesTest(const BitsStorage& data);

/// @brief Binary Matrix Rank Test.
/// @details The focus of the test is the rank of disjoint sub-matrices of the entire sequence. The purpose of this test is
/// to check for linear dependence among fixed length substrings of the original sequence. Note that this test
/// also appears in the DIEHARD battery of tests.
///
/// @param[in] data Class contained sequence for test.
/// @return @ref nistpp::return_t
return_t RankTest(const BitsStorage& data);

/// @brief Discrete Fourier Transform (Spectral) Test.
/// @details The focus of this test is the peak heights in the Discrete Fourier Transform of the sequence. The purpose
/// of this test is to detect periodic features (i.e., repetitive patterns that are near each other) in the tested
/// sequence that would indicate a deviation from the assumption of randomness. The intention is to detect
/// whether the number of peaks exceeding the 95 % threshold is significantly different than 5 %.
///
/// @param[in] data Class contained sequence for test.
/// @return @ref nistpp::return_t
return_t FftTest(const BitsStorage& data);

/// @brief Non-overlapping Template Matching Test.
/// @details The focus of this test is the number of occurrences of pre-specified target strings. The purpose of this
/// test is to detect generators that produce too many occurrences of a given non-periodic (aperiodic) pattern.
/// For this test and for the Overlapping Template Matching test, an m-bit window is used to
/// search for a specific m-bit pattern. If the pattern is not found, the window slides one bit position. If the
/// pattern is found, the window is reset to the bit after the found pattern, and the search resumes.
///
/// @param[in] data Class contained sequence for test.
/// @param[in] m    The length in bits of each template.
/// @param[out] P   Array P-value for all rows in template.
/// @return @ref nistpp::return_t
return_t NonOverlappingTemplateTest(const BitsStorage& data, std::size_t m, std::vector<double>& P);

/// @brief Overlapping Template Matching Test.
/// @details The focus of the Overlapping Template Matching test is the number of occurrences of pre-specified target
/// strings. Both this test and the Non-overlapping Template Matching test of Section 2.7 use an m-bit
/// window to search for a specific m-bit pattern. As with the test in Section 2.7, if the pattern is not found,
/// the window slides one bit position. The difference between this test and the test in Section 2.7 is that
/// when the pattern is found, the window slides only one bit before resuming the search.
///
/// @param[in] data Class contained sequence for test.
/// @param[in] m    The length in bits of each template.
/// @return @ref nistpp::return_t
return_t OverlappingTemplateTest(const BitsStorage& data, std::size_t m);

/// @brief Maurer’s “Universal Statistical”.
/// @details The focus of this test is the number of bits between matching patterns (a measure that is related to the
/// length of a compressed sequence). The purpose of the test is to detect whether or not the sequence can be
/// significantly compressed without loss of information. A significantly compressible sequence is
/// considered to be non-random.
///
/// @param[in] data Class contained sequence for test.
/// @return @ref nistpp::return_t
return_t UniversalTest(const BitsStorage& data);

/// @brief Linear Complexity Test.
/// @details The focus of this test is the length of a linear feedback shift register (LFSR). The purpose of this test is to
/// determine whether or not the sequence is complex enough to be considered random. Random sequences
/// are characterized by longer LFSRs. An LFSR that is too short implies non-randomness.
///
/// @param[in] data Class contained sequence for test.
/// @param[in] M    The length in bits of a block.
/// @return @ref nistpp::return_t
return_t LinearComplexityTest(const BitsStorage& data, std::size_t M);

/// @brief Serial Test.
/// @details The focus of this test is the frequency of all possible overlapping m-bit patterns across the entire
/// sequence. The purpose of this test is to determine whether the number of occurrences of the 2m m-bit
/// overlapping patterns is approximately the same as would be expected for a random sequence. Random
/// sequences have uniformity; that is, every m-bit pattern has the same chance of appearing as every other
/// m-bit pattern. Note that for m = 1, the Serial test is equivalent to the Frequency test.
///
/// @param[in] data Class contained sequence for test.
/// @param[in] M    The length in bits of a block.
/// @param[out] P   Two P-values.
/// @return @ref nistpp::return_t
return_t SerialTest(const BitsStorage &data, std::size_t M, std::array<double, 2> &P);

/// @brief Approximate Entropy Test.
/// @details As with the Serial test, the focus of this test is the frequency of all possible overlapping
/// m-bit patterns across the entire sequence. The purpose of the test is to compare the frequency of
/// overlapping blocks of two consecutive/adjacent lengths (m and m+1) against the expected result for a
/// random sequence.
///
/// @param[in] data Class contained sequence for test.
/// @param[in] M    The length of each block – in this case, the first block length used in the test. m+1 is the
///                 second block length used.
/// @return @ref nistpp::return_t
return_t ApproximateEntropyTest(const BitsStorage& data, std::size_t M);

/// @brief Cumulative Sums (Cusum) Test.
/// @details The focus of this test is the maximal excursion (from zero) of the random walk defined by the cumulative
/// sum of adjusted (-1, +1) digits in the sequence. The purpose of the test is to determine whether the
/// cumulative sum of the partial sequences occurring in the tested sequence is too large or too small relative
/// to the expected behavior of that cumulative sum for random sequences. This cumulative sum may be
/// considered as a random walk. For a random sequence, the excursions of the random walk should be near
/// zero. For certain types of non-random sequences, the excursions of this random walk from zero will be
/// large.
///
/// @param[in] data Class contained sequence for test.
/// @param[out] P   P-values for cusum-forward and cusum-reverse.
/// @return @ref nistpp::return_t
return_t CumulativeSumsTest(const BitsStorage &data, std::array<double, 2> &P);

/// @brief Random Excursions Test.
/// @details The focus of this test is the number of cycles having exactly K visits in a cumulative sum random walk.
/// The cumulative sum random walk is derived from partial sums after the (0,1) sequence is transferred to
/// the appropriate (-1, +1) sequence. A cycle of a random walk consists of a sequence of steps of unit length
/// taken at random that begin at and return to the origin. The purpose of this test is to determine if the
/// number of visits to a particular state within a cycle deviates from what one would expect for a random
/// sequence. This test is actually a series of eight tests (and conclusions), one test and conclusion for each of
/// the states: -4, -3, -2, -1 and +1, +2, +3, +4.
///
/// @param[in]  data Class contained sequence for test.
/// @param[out] P    P-value for 8 state: [-4, -3, -2, -1, 1, 2, 3, 4].
/// @return @ref nistpp::return_t
/// @throw std::runtime_error Input data incorrect for this test.
return_t RandomExcursionsTest(const BitsStorage& data, std::array<double, 8> &P);

/// @brief Random Excursions Variant Test.
/// @details The focus of this test is the total number of times that a particular state is visited (i.e., occurs) in a
/// cumulative sum random walk. The purpose of this test is to detect deviations from the expected number
/// of visits to various states in the random walk. This test is actually a series of eighteen tests (and
/// conclusions), one test and conclusion for each of the states: -9, -8, …, -1 and +1, +2, …, +9.
///
/// @param[in]  data Class contained sequence for test.
/// @param[out] P    P-value for 18 state: [-9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9].
/// @return @ref nistpp::return_t
return_t RandomExcursionsVariantTest(const BitsStorage& data, std::array<double, 18>& P);

} // namespace nistpp

#endif // TESTS_H
