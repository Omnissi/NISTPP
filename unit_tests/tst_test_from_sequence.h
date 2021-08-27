#ifndef TST_TEST_FROM_SEQUENCE_H
#define TST_TEST_FROM_SEQUENCE_H

#include <gtest/gtest.h>

#include <nistpp/bits_storage.h>

#include <nistpp/tests.h>

#include <boost/filesystem.hpp>

#include <cmath>

class nistpp_sequence_test : public ::testing::Test
{
protected:
    void SetUp()
    {
        boost::filesystem::path path(std::string(FILE_PREFIX) + std::string("/data.sha1"));
        if( !boost::filesystem::exists( path ) )
        {
            FAIL() << "File [" << path.c_str() << "] doesn't exist!";
        }
        std::ifstream stream;
        stream.open(path.string(), std::ios::binary);

        if(stream.fail() || stream.bad())
        {
            FAIL() << "Can't open file: " << stream.failbit;
        }

        std::vector<uint8_t> tmp(boost::filesystem::file_size(path));
        std::size_t size = 0;

        while(size < tmp.size())
        {
            auto err = stream.readsome(reinterpret_cast<char*>(tmp.data()), tmp.size());
            if(err < 0)
            {
                FAIL() << "Error read sequence";
            }

            size += err;
        }

        bitsStorage_ = std::make_shared<nistpp::BitsStorage>(tmp);
    }

    void TearDown()
    {}

    bool EqualWithNistPValue(const double& val, const double& PNist)
    {
        return std::fabs(val - PNist) < 1e-6;
    }

    std::shared_ptr<nistpp::BitsStorage>    bitsStorage_;
};

TEST_F(nistpp_sequence_test, frequency)
{
    auto res = nistpp::FrequencyTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.604458));
}

TEST_F(nistpp_sequence_test, block_frequency)
{
    auto res = nistpp::BlockFrequencyTest(*bitsStorage_, 128);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.091517));
}

TEST_F(nistpp_sequence_test, runs)
{
    auto res = nistpp::RunsTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.309757));
}

TEST_F(nistpp_sequence_test, long_runs_of_ones)
{
    auto res = nistpp::LongestRunOfOnesTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.657812));
}

TEST_F(nistpp_sequence_test, rank)
{
    auto res = nistpp::RankTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.577829));
}

TEST_F(nistpp_sequence_test, fft)
{
    auto res = nistpp::FftTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.163062));
}
#endif // TST_TEST_FROM_SEQUENCE_H
