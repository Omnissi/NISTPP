#ifndef TST_LONGEST_RUN_OF_ONES_H
#define TST_LONGEST_RUN_OF_ONES_H

#include <gtest/gtest.h>

#include <nistpp/tests.h>
#include <nistpp/bits_storage.h>
#include <nistpp/sequence_helpers.h>

#include <cmath>

using namespace testing;

TEST(nistpp, longestRunOfOnesTest)
{
    constexpr double passPValue = 0.18060931823971207;
    auto data = nistpp::ConvertStringToSequence("11001100000101010110110001001100111000000000001001"
                                                "00110101010001000100111101011010000000110101111100"
                                                "1100111001101101100010110010");

    nistpp::BitsStorage bits(data);
    auto res = nistpp::LongestRunOfOnesTest(bits);

    EXPECT_TRUE(std::get<0>(res));
    EXPECT_EQ(passPValue, std::get<1>(res));
}

#endif // TST_LONGEST_RUN_OF_ONES_H
