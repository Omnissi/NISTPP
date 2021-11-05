/**
 * @file
 *
 * @author Negodyaev Sergey (negodyaev.sergey@outlook.com)
 * @details This file contains declarations of all tests.
 */
#ifndef TESTS_H
#define TESTS_H

#include <nistpp/bits_storage.h>
#include <nistpp/types.h>

namespace nistpp
{

/// @brief Frequency (Monobit) Test
/// @details The focus of the test is the proportion of zeroes and ones for the entire sequence. The purpose of this test
/// is to determine whether the number of ones and zeros in a sequence are approximately the same as would
/// be expected for a truly random sequence. The test assesses the closeness of the fraction of ones to 1/2, that
/// is, the number of ones and zeroes in a sequence should be about the same. All subsequent tests depend on
/// the passing of this test.
///
/// @param[in] data Class contained sequnce for test.
/// @return @ref nistpp::return_t
return_t FrequencyTest(const BitsStorage& data);

/// @brief Frequency Test within a Block.
/// @details The focus of the test is the proportion of ones within M-bit blocks. The purpose of this test is to determine
/// whether the frequency of ones in an M-bit block is approximately M/2, as would be expected under an
/// assumption of randomness. For block size M=1, this test degenerates to test 1, the Frequency (Monobit) test.
///
/// @param[in] data Class contained sequnce for test.
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
/// @param[in] data Class contained sequnce for test.
/// @return @ref nistpp::return_t
return_t RunsTest(const BitsStorage& data);

/// @brief Test for the Longest Run of Ones in a Block.
/// @details The focus of the test is the longest run of ones within M-bit blocks. The purpose of this test is to
/// determine whether the length of the longest run of ones within the tested sequence is consistent with the
/// length of the longest run of ones that would be expected in a random sequence. Note that an irregularity in
/// the expected length of the longest run of ones implies that there is also an irregularity in the expected
/// length of the longest run of zeroes. Therefore, only a test for ones is necessary.
///
/// @param[in] data Class contained sequnce for test.
/// @return @ref nistpp::return_t
return_t LongestRunOfOnesTest(const BitsStorage& data);

/// @brief Binary Matrix Rank Test.
/// @details The focus of the test is the rank of disjoint sub-matrices of the entire sequence. The purpose of this test is
/// to check for linear dependence among fixed length substrings of the original sequence. Note that this test
/// also appears in the DIEHARD battery of tests.
///
/// @param[in] data Class contained sequnce for test.
/// @return @ref nistpp::return_t
return_t RankTest(const BitsStorage& data);

/// @brief Discrete Fourier Transform (Spectral) Test.
/// @details The focus of this test is the peak heights in the Discrete Fourier Transform of the sequence. The purpose
/// of this test is to detect periodic features (i.e., repetitive patterns that are near each other) in the tested
/// sequence that would indicate a deviation from the assumption of randomness. The intention is to detect
/// whether the number of peaks exceeding the 95 % threshold is significantly different than 5 %.
///
/// @param[in] data Class contained sequnce for test.
/// @return @ref nistpp::return_t
return_t FftTest(const BitsStorage& data);

/// @brief Non-overlapping Template Matching Test.
/// @details The focus of this test is the number of occurrences of pre-specified target strings. The purpose of this
/// test is to detect generators that produce too many occurrences of a given non-periodic (aperiodic) pattern.
/// For this test and for the Overlapping Template Matching test, an m-bit window is used to
/// search for a specific m-bit pattern. If the pattern is not found, the window slides one bit position. If the
/// pattern is found, the window is reset to the bit after the found pattern, and the search resumes.
///
/// @param[in] data Class contained sequnce for test.
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
/// @param[in] data Class contained sequnce for test.
/// @param[in] m    The length in bits of each template.
/// @return @ref nistpp::return_t
return_t OverlappingTemplateTest(const BitsStorage& data, std::size_t m);

} // namespace nistpp

#endif // TESTS_H
