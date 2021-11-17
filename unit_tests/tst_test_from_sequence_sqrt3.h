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
    void SetUp() override
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

    void TearDown() override
    {}

    static bool EqualWithNistPValue(const double& val, const double& PNist)
    {
        return std::fabs(val - PNist) < 1e-6;
    }

    std::shared_ptr<nistpp::BitsStorage>    bitsStorage_{};
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
    EXPECT_TRUE(EqualWithNistPValue(P[1], 0.887225));
    EXPECT_TRUE(EqualWithNistPValue(P[2], 0.328113));
    EXPECT_TRUE(EqualWithNistPValue(P[3], 0.992024));
    EXPECT_TRUE(EqualWithNistPValue(P[4], 0.854501));
    EXPECT_TRUE(EqualWithNistPValue(P[5], 0.273225));
    EXPECT_TRUE(EqualWithNistPValue(P[6], 0.389844));
    EXPECT_TRUE(EqualWithNistPValue(P[7], 0.805707));
    EXPECT_TRUE(EqualWithNistPValue(P[8], 0.238178));
    EXPECT_TRUE(EqualWithNistPValue(P[9], 0.981848));
    EXPECT_TRUE(EqualWithNistPValue(P[10], 0.986376));
    EXPECT_TRUE(EqualWithNistPValue(P[11], 0.688552));
    EXPECT_TRUE(EqualWithNistPValue(P[12], 0.513210));
    EXPECT_TRUE(EqualWithNistPValue(P[13], 0.226998));
    EXPECT_TRUE(EqualWithNistPValue(P[14], 0.753423));
    EXPECT_TRUE(EqualWithNistPValue(P[15], 0.033345));
    EXPECT_TRUE(EqualWithNistPValue(P[16], 0.194836));
    EXPECT_TRUE(EqualWithNistPValue(P[17], 0.253272));
    EXPECT_TRUE(EqualWithNistPValue(P[18], 0.399046));
    EXPECT_TRUE(EqualWithNistPValue(P[19], 0.183303));
    EXPECT_TRUE(EqualWithNistPValue(P[20], 0.891072));
    EXPECT_TRUE(EqualWithNistPValue(P[21], 0.394578));
    EXPECT_TRUE(EqualWithNistPValue(P[22], 0.822593));
    EXPECT_TRUE(EqualWithNistPValue(P[23], 0.946055));
    EXPECT_TRUE(EqualWithNistPValue(P[24], 0.676465));
    EXPECT_TRUE(EqualWithNistPValue(P[25], 0.208579));
    EXPECT_TRUE(EqualWithNistPValue(P[26], 0.549645));
    EXPECT_TRUE(EqualWithNistPValue(P[27], 0.961691));
    EXPECT_TRUE(EqualWithNistPValue(P[28], 0.389116));
    EXPECT_TRUE(EqualWithNistPValue(P[29], 0.362742));
    EXPECT_TRUE(EqualWithNistPValue(P[30], 0.684970));
    EXPECT_TRUE(EqualWithNistPValue(P[31], 0.199787));
    EXPECT_TRUE(EqualWithNistPValue(P[32], 0.723818));
    EXPECT_TRUE(EqualWithNistPValue(P[33], 0.824072));
    EXPECT_TRUE(EqualWithNistPValue(P[34], 0.703254));
    EXPECT_TRUE(EqualWithNistPValue(P[35], 0.888162));
    EXPECT_TRUE(EqualWithNistPValue(P[36], 0.321256));
    EXPECT_TRUE(EqualWithNistPValue(P[37], 0.122400));
    EXPECT_TRUE(EqualWithNistPValue(P[38], 0.924783));
    EXPECT_TRUE(EqualWithNistPValue(P[39], 0.064828));
    EXPECT_TRUE(EqualWithNistPValue(P[40], 0.506507));
    EXPECT_TRUE(EqualWithNistPValue(P[41], 0.377596));
    EXPECT_TRUE(EqualWithNistPValue(P[42], 0.713390));
    EXPECT_TRUE(EqualWithNistPValue(P[43], 0.740004));
    EXPECT_TRUE(EqualWithNistPValue(P[44], 0.497652));
    EXPECT_TRUE(EqualWithNistPValue(P[45], 0.256037));
    EXPECT_TRUE(EqualWithNistPValue(P[46], 0.466285));
    EXPECT_TRUE(EqualWithNistPValue(P[47], 0.202168));
    EXPECT_TRUE(EqualWithNistPValue(P[48], 0.227588));
    EXPECT_TRUE(EqualWithNistPValue(P[49], 0.131734));
    EXPECT_TRUE(EqualWithNistPValue(P[50], 0.158560));
    EXPECT_TRUE(EqualWithNistPValue(P[51], 0.109698));
    EXPECT_TRUE(EqualWithNistPValue(P[52], 0.574424));
    EXPECT_TRUE(EqualWithNistPValue(P[53], 0.460508));
    EXPECT_TRUE(EqualWithNistPValue(P[54], 0.436346));
    EXPECT_TRUE(EqualWithNistPValue(P[55], 0.816963));
    EXPECT_TRUE(EqualWithNistPValue(P[56], 0.417371));
    EXPECT_TRUE(EqualWithNistPValue(P[57], 0.781368));
    EXPECT_TRUE(EqualWithNistPValue(P[58], 0.957793));
    EXPECT_TRUE(EqualWithNistPValue(P[59], 0.747372));
    EXPECT_TRUE(EqualWithNistPValue(P[60], 0.402440));
    EXPECT_TRUE(EqualWithNistPValue(P[61], 0.098035));
    EXPECT_TRUE(EqualWithNistPValue(P[62], 0.767256));
    EXPECT_TRUE(EqualWithNistPValue(P[63], 0.485520));
    EXPECT_TRUE(EqualWithNistPValue(P[64], 0.666936));
    EXPECT_TRUE(EqualWithNistPValue(P[65], 0.354678));
    EXPECT_TRUE(EqualWithNistPValue(P[66], 0.596152));
    EXPECT_TRUE(EqualWithNistPValue(P[67], 0.134610));
    EXPECT_TRUE(EqualWithNistPValue(P[68], 0.680591));
    EXPECT_TRUE(EqualWithNistPValue(P[69], 0.939800));
    EXPECT_TRUE(EqualWithNistPValue(P[70], 0.901293));
    EXPECT_TRUE(EqualWithNistPValue(P[71], 0.079350));
    EXPECT_TRUE(EqualWithNistPValue(P[72], 0.176848));
    EXPECT_TRUE(EqualWithNistPValue(P[73], 0.071205));
    EXPECT_TRUE(EqualWithNistPValue(P[74], 0.424530));
    EXPECT_TRUE(EqualWithNistPValue(P[75], 0.247669));
    EXPECT_TRUE(EqualWithNistPValue(P[76], 0.273731));
    EXPECT_TRUE(EqualWithNistPValue(P[77], 0.788742));
    EXPECT_TRUE(EqualWithNistPValue(P[78], 0.598447));
    EXPECT_TRUE(EqualWithNistPValue(P[79], 0.004249));
    EXPECT_TRUE(EqualWithNistPValue(P[80], 0.322169));
    EXPECT_TRUE(EqualWithNistPValue(P[81], 0.570237));
    EXPECT_TRUE(EqualWithNistPValue(P[82], 0.308038));
    EXPECT_TRUE(EqualWithNistPValue(P[83], 0.957103));
    EXPECT_TRUE(EqualWithNistPValue(P[84], 0.930850));
    EXPECT_TRUE(EqualWithNistPValue(P[85], 0.198524));
    EXPECT_TRUE(EqualWithNistPValue(P[86], 0.041746));
    EXPECT_TRUE(EqualWithNistPValue(P[87], 0.555491));
    EXPECT_TRUE(EqualWithNistPValue(P[88], 0.839062));
    EXPECT_TRUE(EqualWithNistPValue(P[89], 0.137494));
    EXPECT_TRUE(EqualWithNistPValue(P[90], 0.421857));
    EXPECT_TRUE(EqualWithNistPValue(P[91], 0.849809));
    EXPECT_TRUE(EqualWithNistPValue(P[92], 0.135209));
    EXPECT_TRUE(EqualWithNistPValue(P[93], 0.768058));
    EXPECT_TRUE(EqualWithNistPValue(P[94], 0.004231));
    EXPECT_TRUE(EqualWithNistPValue(P[95], 0.605397));
    EXPECT_TRUE(EqualWithNistPValue(P[96], 0.161513));
    EXPECT_TRUE(EqualWithNistPValue(P[97], 0.028669));
    EXPECT_TRUE(EqualWithNistPValue(P[98], 0.785612));
    EXPECT_TRUE(EqualWithNistPValue(P[99], 0.319326));
    EXPECT_TRUE(EqualWithNistPValue(P[100], 0.917332));
    EXPECT_TRUE(EqualWithNistPValue(P[101], 0.289763));
    EXPECT_TRUE(EqualWithNistPValue(P[102], 0.402166));
    EXPECT_TRUE(EqualWithNistPValue(P[103], 0.977245));
    EXPECT_TRUE(EqualWithNistPValue(P[104], 0.547568));
    EXPECT_TRUE(EqualWithNistPValue(P[105], 0.666999));
    EXPECT_TRUE(EqualWithNistPValue(P[106], 0.658943));
    EXPECT_TRUE(EqualWithNistPValue(P[107], 0.310540));
    EXPECT_TRUE(EqualWithNistPValue(P[108], 0.616125));
    EXPECT_TRUE(EqualWithNistPValue(P[109], 0.914996));
    EXPECT_TRUE(EqualWithNistPValue(P[110], 0.875309));
    EXPECT_TRUE(EqualWithNistPValue(P[111], 0.272204));
    EXPECT_TRUE(EqualWithNistPValue(P[112], 0.188305));
    EXPECT_TRUE(EqualWithNistPValue(P[113], 0.146763));
    EXPECT_TRUE(EqualWithNistPValue(P[114], 0.055417));
    EXPECT_TRUE(EqualWithNistPValue(P[115], 0.256215));
    EXPECT_TRUE(EqualWithNistPValue(P[116], 0.969386));
    EXPECT_TRUE(EqualWithNistPValue(P[117], 0.009763));
    EXPECT_TRUE(EqualWithNistPValue(P[118], 0.079618));
    EXPECT_TRUE(EqualWithNistPValue(P[119], 0.270477));
    EXPECT_TRUE(EqualWithNistPValue(P[120], 0.315155));
    EXPECT_TRUE(EqualWithNistPValue(P[121], 0.313453));
    EXPECT_TRUE(EqualWithNistPValue(P[122], 0.906179));
    EXPECT_TRUE(EqualWithNistPValue(P[123], 0.827663));
    EXPECT_TRUE(EqualWithNistPValue(P[124], 0.449617));
    EXPECT_TRUE(EqualWithNistPValue(P[125], 0.507883));
    EXPECT_TRUE(EqualWithNistPValue(P[126], 0.218126));
    EXPECT_TRUE(EqualWithNistPValue(P[127], 0.766346));
    EXPECT_TRUE(EqualWithNistPValue(P[128], 0.919958));
    EXPECT_TRUE(EqualWithNistPValue(P[129], 0.126258));
    EXPECT_TRUE(EqualWithNistPValue(P[130], 0.096635));
    EXPECT_TRUE(EqualWithNistPValue(P[131], 0.499688));
    EXPECT_TRUE(EqualWithNistPValue(P[132], 0.181426));
    EXPECT_TRUE(EqualWithNistPValue(P[133], 0.001628));
    EXPECT_TRUE(EqualWithNistPValue(P[134], 0.519370));
    EXPECT_TRUE(EqualWithNistPValue(P[135], 0.203758));
    EXPECT_TRUE(EqualWithNistPValue(P[136], 0.007882));
    EXPECT_TRUE(EqualWithNistPValue(P[137], 0.947168));
    EXPECT_TRUE(EqualWithNistPValue(P[138], 0.759773));
    EXPECT_TRUE(EqualWithNistPValue(P[139], 0.424877));
    EXPECT_TRUE(EqualWithNistPValue(P[140], 0.101836));
    EXPECT_TRUE(EqualWithNistPValue(P[141], 0.149678));
    EXPECT_TRUE(EqualWithNistPValue(P[142], 0.267217));
    EXPECT_TRUE(EqualWithNistPValue(P[143], 0.754730));
    EXPECT_TRUE(EqualWithNistPValue(P[144], 0.168805));
    EXPECT_TRUE(EqualWithNistPValue(P[145], 0.600663));
    EXPECT_TRUE(EqualWithNistPValue(P[146], 0.248307));
    EXPECT_TRUE(EqualWithNistPValue(P[147], 0.071205));
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
    std::array<double, 2> P;
    auto res = nistpp::SerialTest(*bitsStorage_, 16, P);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.180826));
}

TEST_F(nistpp_sequence_test_sqrt3, Approximate)
{
    auto res = nistpp::ApproximateEntropyTest(*bitsStorage_, 10);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.154929));
}

TEST_F(nistpp_sequence_test_sqrt3, Cusum)
{
    std::array<double, 2> P;
    auto res = nistpp::CumulativeSumsTest(*bitsStorage_, P);

    EXPECT_TRUE(EqualWithNistPValue(std::get<1>(res), 0.605167));
}

TEST_F(nistpp_sequence_test_sqrt3, RandomExcursions)
{
    std::array<double, 8> P{};
    nistpp::RandomExcursionsTest(*bitsStorage_, P);

    EXPECT_TRUE(EqualWithNistPValue(P[4], 0.783283));
}

TEST_F(nistpp_sequence_test_sqrt3, RandomExcursionsVariant)
{
    std::array<double, 18> P{};
    nistpp::RandomExcursionsVariantTest(*bitsStorage_, P);

    EXPECT_TRUE(EqualWithNistPValue(P[8], 0.155066));
}


#endif //NISTPP_TST_TEST_FROM_SEQUENCE_SQRT3_H
