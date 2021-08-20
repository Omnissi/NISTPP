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

} // namespace nistpp

#endif // TESTS_H
