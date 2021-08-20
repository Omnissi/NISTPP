#ifndef TST_RUNS_H
#define TST_RUNS_H

#include <gtest/gtest.h>

#include <nistpp/tests.h>
#include <nistpp/bits_storage.h>
#include <nistpp/sequence_helpers.h>

#include <cmath>

using namespace testing;

TEST(nistpp, runs)
{
    constexpr double passPValue = 0.6953317934158357;
    auto data = nistpp::ConvertStringToSequence("11001001000011111101101010100010001000010110100011"
                                                "000010001101001100010011000110011000101000101110000000");

    nistpp::BitsStorage bits(data);
    auto res = nistpp::RunsTest(bits);

    EXPECT_TRUE(std::get<0>(res));
    EXPECT_EQ(passPValue, std::get<1>(res));
}

#endif // TST_RUNS_H
