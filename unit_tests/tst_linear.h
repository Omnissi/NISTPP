#ifndef TST_LINEAR_H
#define TST_LINEAR_H

#include <gtest/gtest.h>

#include <nistpp/tests.h>
#include <nistpp/bits_storage.h>
#include <nistpp/sequence_helpers.h>

#include <boost/filesystem.hpp>

#include <cmath>
#include <string>

using namespace testing;

TEST(nistpp, serial)
{
    boost::filesystem::path path(std::string(FILE_PREFIX) + std::string("/data.e"));
    if(!boost::filesystem::exists(path))
    {
        FAIL() << "File [" << path.c_str() << "] doesn't exist!";
    }
    std::ifstream stream;
    stream.open(path.string());

    if(stream.fail() || stream.bad())
    {
        FAIL() << "Can't open file: " << stream.failbit;
    }

    std::string res;
    std::string tmp;

    while(std::getline(stream, tmp) && static_cast<double>(res.size()) <= 1e6)
    {
        res += tmp;
    }

    res.resize(1e6);
    auto data = nistpp::ConvertStringToSequence(res);

    nistpp::BitsStorage bits(data);

    std::array<double, 2> P{};
    nistpp::SerialTest(bits, 2, P);

    EXPECT_TRUE(std::fabs(P[0] - 0.843764) < 1e-6);
    EXPECT_TRUE(std::fabs(P[1] - 0.561915) < 1e-6);
}

#endif // TST_LINEAR_H
