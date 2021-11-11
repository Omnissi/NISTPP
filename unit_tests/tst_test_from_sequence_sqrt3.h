#ifndef NISTPP_TST_TEST_FROM_SEQUENCE_SQRT3_H
#define NISTPP_TST_TEST_FROM_SEQUENCE_SQRT3_H

#include <gtest/gtest.h>

#include <nistpp/bits_storage.h>

#include <nistpp/tests.h>

#include <boost/filesystem.hpp>

#include <cmath>
#include <nistpp/sequence_helpers.h>

class nistpp_sequence_test_sqrt3 : public ::testing::Test
{
protected:
    void SetUp()
    {
        boost::filesystem::path path(std::string(FILE_PREFIX) + std::string("/data.sqrt3"));
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

        while(std::getline(stream, tmp))
        {
            res += tmp;
        }

        auto data = nistpp::ConvertStringToSequence(res);

        bitsStorage_ = std::make_shared<nistpp::BitsStorage>(data);
    }

    void TearDown()
    {}

    bool EqualWithNistPValue(const double& val, const double& PNist)
    {
        return std::fabs(val - PNist) < 1e-6;
    }

    std::shared_ptr<nistpp::BitsStorage>    bitsStorage_;
};

TEST_F(nistpp_sequence_test_sqrt3, frequency)
{
    auto res = nistpp::FrequencyTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.546820));
}

TEST_F(nistpp_sequence_test_sqrt3, block_frequency)
{
    auto res = nistpp::BlockFrequencyTest(*bitsStorage_, 128);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.473925));
}

TEST_F(nistpp_sequence_test_sqrt3, runs)
{
    auto res = nistpp::RunsTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.229863));
}

TEST_F(nistpp_sequence_test_sqrt3, long_runs_of_ones)
{
    auto res = nistpp::LongestRunOfOnesTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.446726));
}

TEST_F(nistpp_sequence_test_sqrt3, rank)
{
    auto res = nistpp::RankTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.348786));
}

TEST_F(nistpp_sequence_test_sqrt3, fft)
{
    auto res = nistpp::FftTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.509824));
}

TEST_F(nistpp_sequence_test_sqrt3, nonOverlapping)
{
    std::vector<double> P;
    nistpp::NonOverlappingTemplateTest(*bitsStorage_, 9, P);

    EXPECT_TRUE(EqualWithNistPValue(P[0], 0.424530));
}

TEST_F(nistpp_sequence_test_sqrt3, Overlapping)
{
    auto res = nistpp::OverlappingTemplateTest(*bitsStorage_, 9);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.070981));
}

TEST_F(nistpp_sequence_test_sqrt3, Universal)
{
    auto res = nistpp::UniversalTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.150578));
}

TEST_F(nistpp_sequence_test_sqrt3, Linear)
{
    auto res = nistpp::LinearComplexityTest(*bitsStorage_, 500);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.341994));
}

TEST_F(nistpp_sequence_test_sqrt3, Serial)
{
    auto res = nistpp::SerialTest(*bitsStorage_, 16);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.180826));
}

TEST_F(nistpp_sequence_test_sqrt3, Approximate)
{
    auto res = nistpp::ApproximateEntropyTest(*bitsStorage_, 10);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.154929));
}

TEST_F(nistpp_sequence_test_sqrt3, Cusum)
{
    auto res = nistpp::CumulativeSumsTest(*bitsStorage_);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.605167));
}

TEST_F(nistpp_sequence_test_sqrt3, RandomExcursions)
{
    std::array<double, 8> P;
    nistpp::RandomExcursionsTest(*bitsStorage_, P);

    EXPECT_TRUE(EqualWithNistPValue(P[4], 0.783283));
}

TEST_F(nistpp_sequence_test_sqrt3, RandomExcursionsVariant)
{
    std::array<double, 18> P;
    nistpp::RandomExcursionsVariantTest(*bitsStorage_, P);

    EXPECT_TRUE(EqualWithNistPValue(P[8], 0.155066));
}


#endif //NISTPP_TST_TEST_FROM_SEQUENCE_SQRT3_H
