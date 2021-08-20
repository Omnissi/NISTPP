#ifndef TST_FREQUENCY_H
#define TST_FREQUENCY_H

#include <gtest/gtest.h>

#include <nistpp/tests.h>
#include <nistpp/bits_storage.h>
#include <nistpp/sequence_helpers.h>

#include <cmath>

using namespace testing;

TEST(nistpp, freuqency)
{
    constexpr double passPValue = 0.11666446478102338;
    auto data = nistpp::ConvertStringToSequence("1100100100001111110110101010001000100001011010001100"
                                                "0010001101001100010011000110011000101000101110001100");

    nistpp::BitsStorage bits(data);
    auto res = nistpp::FrequencyTest(bits);

    EXPECT_TRUE(std::get<0>(res));
    EXPECT_EQ(passPValue, std::get<1>(res));
}

#endif // TST_FREQUENCY_H
