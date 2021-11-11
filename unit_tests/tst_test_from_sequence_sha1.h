#ifndef TST_TEST_FROM_SEQUENCE_H
#define TST_TEST_FROM_SEQUENCE_H

#include <gtest/gtest.h>

#include <nistpp/bits_storage.h>

#include <nistpp/tests.h>

#include <boost/filesystem.hpp>

#include <cmath>

class nistpp_sequence_test_sha1 : public ::testing::Test
{
protected:
    void SetUp()
    {
        boost::filesystem::path path(std::string(FILE_PREFIX) + std::string("/data.sha1"));
        if(!boost::filesystem::exists(path))
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
            auto err = stream.readsome(reinterpret_cast<char*>(tmp.data()) + size, tmp.size() - size);
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

TEST_F(nistpp_sequence_test_sha1, frequency)
{
    auto res = nistpp::FrequencyTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.604458));
}

TEST_F(nistpp_sequence_test_sha1, block_frequency)
{
    auto res = nistpp::BlockFrequencyTest(*bitsStorage_, 128);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.091517));
}

TEST_F(nistpp_sequence_test_sha1, runs)
{
    auto res = nistpp::RunsTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.309757));
}

TEST_F(nistpp_sequence_test_sha1, long_runs_of_ones)
{
    auto res = nistpp::LongestRunOfOnesTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.657812));
}

TEST_F(nistpp_sequence_test_sha1, rank)
{
    auto res = nistpp::RankTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.577829));
}

TEST_F(nistpp_sequence_test_sha1, fft)
{
    auto res = nistpp::FftTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.163062));
}

TEST_F(nistpp_sequence_test_sha1, nonOverlapping)
{
    std::vector<double> P;
    auto res = nistpp::NonOverlappingTemplateTest(*bitsStorage_, 9, P);

    EXPECT_TRUE(EqualWithNistPValue(P[0], 0.496601));
}

TEST_F(nistpp_sequence_test_sha1, Overlapping)
{
    auto res = nistpp::OverlappingTemplateTest(*bitsStorage_, 9);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.339426));
}

TEST_F(nistpp_sequence_test_sha1, Universal)
{
    auto res = nistpp::UniversalTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.411079));
}

TEST_F(nistpp_sequence_test_sha1, Linear)
{
    auto res = nistpp::LinearComplexityTest(*bitsStorage_, 500);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.309412));
}

TEST_F(nistpp_sequence_test_sha1, Serial)
{
    auto res = nistpp::SerialTest(*bitsStorage_, 16);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.760793));
}

TEST_F(nistpp_sequence_test_sha1, Approximate)
{
    auto res = nistpp::ApproximateEntropyTest(*bitsStorage_, 10);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.982885));
}

TEST_F(nistpp_sequence_test_sha1, Cusum)
{
    auto res = nistpp::CumulativeSumsTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.451231));
}

TEST_F(nistpp_sequence_test_sha1, RandomExcursions)
{
    std::array<double, 8> P;
    EXPECT_THROW(nistpp::RandomExcursionsTest(*bitsStorage_, P), std::runtime_error);
}

TEST_F(nistpp_sequence_test_sha1, RandomExcursionsVariant)
{
    std::array<double, 18> P;
    EXPECT_THROW(nistpp::RandomExcursionsVariantTest(*bitsStorage_, P), std::runtime_error);
}

#endif // TST_TEST_FROM_SEQUENCE_H
