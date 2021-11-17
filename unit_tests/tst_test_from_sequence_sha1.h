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
    void SetUp() override
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

    void TearDown() override
    {}

    static bool EqualWithNistPValue(const double& val, const double& PNist)
    {
        return std::fabs(val - PNist) < 1e-6;
    }

    std::shared_ptr<nistpp::BitsStorage>    bitsStorage_{};
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
    nistpp::NonOverlappingTemplateTest(*bitsStorage_, 9, P);

    EXPECT_TRUE(EqualWithNistPValue(P[0], 0.496601));
    EXPECT_TRUE(EqualWithNistPValue(P[1], 0.421114));
    EXPECT_TRUE(EqualWithNistPValue(P[2], 0.762228));
    EXPECT_TRUE(EqualWithNistPValue(P[3], 0.313147));
    EXPECT_TRUE(EqualWithNistPValue(P[4], 0.267553));
    EXPECT_TRUE(EqualWithNistPValue(P[5], 0.569810));
    EXPECT_TRUE(EqualWithNistPValue(P[6], 0.918072));
    EXPECT_TRUE(EqualWithNistPValue(P[7], 0.314578));
    EXPECT_TRUE(EqualWithNistPValue(P[8], 0.987461));
    EXPECT_TRUE(EqualWithNistPValue(P[9], 0.129969));
    EXPECT_TRUE(EqualWithNistPValue(P[10], 0.001239));
    EXPECT_TRUE(EqualWithNistPValue(P[11], 0.982678));
    EXPECT_TRUE(EqualWithNistPValue(P[12], 0.296665));
    EXPECT_TRUE(EqualWithNistPValue(P[13], 0.502599));
    EXPECT_TRUE(EqualWithNistPValue(P[14], 0.897621));
    EXPECT_TRUE(EqualWithNistPValue(P[15], 0.688362));
    EXPECT_TRUE(EqualWithNistPValue(P[16], 0.943372));
    EXPECT_TRUE(EqualWithNistPValue(P[17], 0.178395));
    EXPECT_TRUE(EqualWithNistPValue(P[18], 0.129886));
    EXPECT_TRUE(EqualWithNistPValue(P[19], 0.679505));
    EXPECT_TRUE(EqualWithNistPValue(P[20], 0.181677));
    EXPECT_TRUE(EqualWithNistPValue(P[21], 0.151164));
    EXPECT_TRUE(EqualWithNistPValue(P[22], 0.045834));
    EXPECT_TRUE(EqualWithNistPValue(P[23], 0.495715));
    EXPECT_TRUE(EqualWithNistPValue(P[24], 0.534630));
    EXPECT_TRUE(EqualWithNistPValue(P[25], 0.173990));
    EXPECT_TRUE(EqualWithNistPValue(P[26], 0.808597));
    EXPECT_TRUE(EqualWithNistPValue(P[27], 0.431466));
    EXPECT_TRUE(EqualWithNistPValue(P[28], 0.021651));
    EXPECT_TRUE(EqualWithNistPValue(P[29], 0.309800));
    EXPECT_TRUE(EqualWithNistPValue(P[30], 0.869885));
    EXPECT_TRUE(EqualWithNistPValue(P[31], 0.839431));
    EXPECT_TRUE(EqualWithNistPValue(P[32], 0.601457));
    EXPECT_TRUE(EqualWithNistPValue(P[33], 0.484811));
    EXPECT_TRUE(EqualWithNistPValue(P[34], 0.837127));
    EXPECT_TRUE(EqualWithNistPValue(P[35], 0.636475));
    EXPECT_TRUE(EqualWithNistPValue(P[36], 0.178341));
    EXPECT_TRUE(EqualWithNistPValue(P[37], 0.076941));
    EXPECT_TRUE(EqualWithNistPValue(P[38], 0.123272));
    EXPECT_TRUE(EqualWithNistPValue(P[39], 0.224647));
    EXPECT_TRUE(EqualWithNistPValue(P[40], 0.793566));
    EXPECT_TRUE(EqualWithNistPValue(P[41], 0.006495));
    EXPECT_TRUE(EqualWithNistPValue(P[42], 0.081035));
    EXPECT_TRUE(EqualWithNistPValue(P[43], 0.496490));
    EXPECT_TRUE(EqualWithNistPValue(P[44], 0.166833));
    EXPECT_TRUE(EqualWithNistPValue(P[45], 0.049395));
    EXPECT_TRUE(EqualWithNistPValue(P[46], 0.760311));
    EXPECT_TRUE(EqualWithNistPValue(P[47], 0.807429));
    EXPECT_TRUE(EqualWithNistPValue(P[48], 0.368466));
    EXPECT_TRUE(EqualWithNistPValue(P[49], 0.219526));
    EXPECT_TRUE(EqualWithNistPValue(P[50], 0.392200));
    EXPECT_TRUE(EqualWithNistPValue(P[51], 0.348746));
    EXPECT_TRUE(EqualWithNistPValue(P[52], 0.911163));
    EXPECT_TRUE(EqualWithNistPValue(P[53], 0.174523));
    EXPECT_TRUE(EqualWithNistPValue(P[54], 0.868871));
    EXPECT_TRUE(EqualWithNistPValue(P[55], 0.792262));
    EXPECT_TRUE(EqualWithNistPValue(P[56], 0.338091));
    EXPECT_TRUE(EqualWithNistPValue(P[57], 0.448151));
    EXPECT_TRUE(EqualWithNistPValue(P[58], 0.883285));
    EXPECT_TRUE(EqualWithNistPValue(P[59], 0.969115));
    EXPECT_TRUE(EqualWithNistPValue(P[60], 0.311638));
    EXPECT_TRUE(EqualWithNistPValue(P[61], 0.504716));
    EXPECT_TRUE(EqualWithNistPValue(P[62], 0.776657));
    EXPECT_TRUE(EqualWithNistPValue(P[63], 0.326533));
    EXPECT_TRUE(EqualWithNistPValue(P[64], 0.041505));
    EXPECT_TRUE(EqualWithNistPValue(P[65], 0.463134));
    EXPECT_TRUE(EqualWithNistPValue(P[66], 0.262498));
    EXPECT_TRUE(EqualWithNistPValue(P[67], 0.108090));
    EXPECT_TRUE(EqualWithNistPValue(P[68], 0.756011));
    EXPECT_TRUE(EqualWithNistPValue(P[69], 0.673945));
    EXPECT_TRUE(EqualWithNistPValue(P[70], 0.948873));
    EXPECT_TRUE(EqualWithNistPValue(P[71], 0.830551));
    EXPECT_TRUE(EqualWithNistPValue(P[72], 0.958124));
    EXPECT_TRUE(EqualWithNistPValue(P[73], 0.541378));
    EXPECT_TRUE(EqualWithNistPValue(P[74], 0.496601));
    EXPECT_TRUE(EqualWithNistPValue(P[75], 0.867577));
    EXPECT_TRUE(EqualWithNistPValue(P[76], 0.896375));
    EXPECT_TRUE(EqualWithNistPValue(P[77], 0.773545));
    EXPECT_TRUE(EqualWithNistPValue(P[78], 0.679505));
    EXPECT_TRUE(EqualWithNistPValue(P[79], 0.610306));
    EXPECT_TRUE(EqualWithNistPValue(P[80], 0.097604));
    EXPECT_TRUE(EqualWithNistPValue(P[81], 0.964100));
    EXPECT_TRUE(EqualWithNistPValue(P[82], 0.941699));
    EXPECT_TRUE(EqualWithNistPValue(P[83], 0.540118));
    EXPECT_TRUE(EqualWithNistPValue(P[84], 0.651542));
    EXPECT_TRUE(EqualWithNistPValue(P[85], 0.018327));
    EXPECT_TRUE(EqualWithNistPValue(P[86], 0.918295));
    EXPECT_TRUE(EqualWithNistPValue(P[87], 0.932746));
    EXPECT_TRUE(EqualWithNistPValue(P[88], 0.574240));
    EXPECT_TRUE(EqualWithNistPValue(P[89], 0.013149));
    EXPECT_TRUE(EqualWithNistPValue(P[90], 0.726873));
    EXPECT_TRUE(EqualWithNistPValue(P[91], 0.510420));
    EXPECT_TRUE(EqualWithNistPValue(P[92], 0.654152));
    EXPECT_TRUE(EqualWithNistPValue(P[93], 0.012088));
    EXPECT_TRUE(EqualWithNistPValue(P[94], 0.284879));
    EXPECT_TRUE(EqualWithNistPValue(P[95], 0.297717));
    EXPECT_TRUE(EqualWithNistPValue(P[96], 0.402032));
    EXPECT_TRUE(EqualWithNistPValue(P[97], 0.731965));
    EXPECT_TRUE(EqualWithNistPValue(P[98], 0.176450));
    EXPECT_TRUE(EqualWithNistPValue(P[99], 0.049447));
    EXPECT_TRUE(EqualWithNistPValue(P[100], 0.965418));
    EXPECT_TRUE(EqualWithNistPValue(P[101], 0.850475));
    EXPECT_TRUE(EqualWithNistPValue(P[102], 0.008562));
    EXPECT_TRUE(EqualWithNistPValue(P[103], 0.062257));
    EXPECT_TRUE(EqualWithNistPValue(P[104], 0.076057));
    EXPECT_TRUE(EqualWithNistPValue(P[105], 0.464314));
    EXPECT_TRUE(EqualWithNistPValue(P[106], 0.095872));
    EXPECT_TRUE(EqualWithNistPValue(P[107], 0.896375));
    EXPECT_TRUE(EqualWithNistPValue(P[108], 0.684585));
    EXPECT_TRUE(EqualWithNistPValue(P[109], 0.616808));
    EXPECT_TRUE(EqualWithNistPValue(P[110], 0.050685));
    EXPECT_TRUE(EqualWithNistPValue(P[111], 0.821608));
    EXPECT_TRUE(EqualWithNistPValue(P[112], 0.462706));
    EXPECT_TRUE(EqualWithNistPValue(P[113], 0.157657));
    EXPECT_TRUE(EqualWithNistPValue(P[114], 0.046889));
    EXPECT_TRUE(EqualWithNistPValue(P[115], 0.825328));
    EXPECT_TRUE(EqualWithNistPValue(P[116], 0.143763));
    EXPECT_TRUE(EqualWithNistPValue(P[117], 0.738651));
    EXPECT_TRUE(EqualWithNistPValue(P[118], 0.755218));
    EXPECT_TRUE(EqualWithNistPValue(P[119], 0.343883));
    EXPECT_TRUE(EqualWithNistPValue(P[120], 0.713373));
    EXPECT_TRUE(EqualWithNistPValue(P[121], 0.169520));
    EXPECT_TRUE(EqualWithNistPValue(P[122], 0.954230));
    EXPECT_TRUE(EqualWithNistPValue(P[123], 0.692841));
    EXPECT_TRUE(EqualWithNistPValue(P[124], 0.424685));
    EXPECT_TRUE(EqualWithNistPValue(P[125], 0.664824));
    EXPECT_TRUE(EqualWithNistPValue(P[126], 0.119901));
    EXPECT_TRUE(EqualWithNistPValue(P[127], 0.457472));
    EXPECT_TRUE(EqualWithNistPValue(P[128], 0.159428));
    EXPECT_TRUE(EqualWithNistPValue(P[129], 0.683404));
    EXPECT_TRUE(EqualWithNistPValue(P[130], 0.531666));
    EXPECT_TRUE(EqualWithNistPValue(P[131], 0.824297));
    EXPECT_TRUE(EqualWithNistPValue(P[132], 0.108906));
    EXPECT_TRUE(EqualWithNistPValue(P[133], 0.725597));
    EXPECT_TRUE(EqualWithNistPValue(P[134], 0.271241));
    EXPECT_TRUE(EqualWithNistPValue(P[135], 0.714774));
    EXPECT_TRUE(EqualWithNistPValue(P[136], 0.783178));
    EXPECT_TRUE(EqualWithNistPValue(P[137], 0.498486));
    EXPECT_TRUE(EqualWithNistPValue(P[138], 0.924418));
    EXPECT_TRUE(EqualWithNistPValue(P[139], 0.575407));
    EXPECT_TRUE(EqualWithNistPValue(P[140], 0.494166));
    EXPECT_TRUE(EqualWithNistPValue(P[141], 0.967612));
    EXPECT_TRUE(EqualWithNistPValue(P[142], 0.458325));
    EXPECT_TRUE(EqualWithNistPValue(P[143], 0.927039));
    EXPECT_TRUE(EqualWithNistPValue(P[144], 0.105639));
    EXPECT_TRUE(EqualWithNistPValue(P[145], 0.017200));
    EXPECT_TRUE(EqualWithNistPValue(P[146], 0.390641));
    EXPECT_TRUE(EqualWithNistPValue(P[147], 0.541378));
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
    std::array<double, 2> P;
    auto res = nistpp::SerialTest(*bitsStorage_, 16, P);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.760793));
}

TEST_F(nistpp_sequence_test_sha1, Approximate)
{
    auto res = nistpp::ApproximateEntropyTest(*bitsStorage_, 10);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.982885));
}

TEST_F(nistpp_sequence_test_sha1, Cusum)
{
    std::array<double, 2> P;
    nistpp::CumulativeSumsTest(*bitsStorage_, P);

    EXPECT_TRUE(EqualWithNistPValue(P[0], 0.451231));
    EXPECT_TRUE(EqualWithNistPValue(P[1], 0.550134));
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
