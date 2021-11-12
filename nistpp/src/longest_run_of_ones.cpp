#include <nistpp/tests.h>

#include <cmath>
#include <numeric>
#include <stdexcept>

#include "help_function.h"

#include <boost/math/special_functions/gamma.hpp>

namespace nistpp
{

constexpr double  threshold     = 0.01;
constexpr size_t  minimalBits   = 128;

return_t LongestRunOfOnesTest(const BitsStorage& data)
{
    const auto numberOfBit = data.NumberOfBits();

    if(numberOfBit < minimalBits)
    {
        throw std::logic_error(std::to_string(minimalBits) + " bit required! In storage: " + std::to_string(numberOfBit));
    }

    size_t K;
    size_t M;
    std::array<size_t, 7> V;
    std::array<size_t, 7> nu = {};
    std::array<double, 7> pi;

    if(numberOfBit < 6272)
    {
        K = 3;
        M = 8;
        V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
        pi[0] = 0.21484375;
        pi[1] = 0.3671875;
        pi[2] = 0.23046875;
        pi[3] = 0.1875;
    }
    else if(numberOfBit < 750000)
    {
        K = 5;
        M = 128;
        V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
        pi[0] = 0.1174035788;
        pi[1] = 0.242955959;
        pi[2] = 0.249363483;
        pi[3] = 0.17517706;
        pi[4] = 0.102701071;
        pi[5] = 0.112398847;
    }
    else
    {
        K = 6;
        M = 10000;
        V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
        pi[0] = 0.0882;
        pi[1] = 0.2092;
        pi[2] = 0.2483;
        pi[3] = 0.1933;
        pi[4] = 0.1208;
        pi[5] = 0.0675;
        pi[6] = 0.0727;
    }

    const size_t numberOfBlock = numberOfBit/M;
    for(size_t i = 0; i < numberOfBlock; ++i)
    {
        size_t maxRuns      = 0;
        size_t currentRuns  = 0;

        for(size_t j = 0; j < M; ++j)
        {
            if(data[i*M +j])
            {
                ++currentRuns;
                if(currentRuns > maxRuns)
                {
                    maxRuns = currentRuns;
                }
            }
            else
            {
                currentRuns = 0;
            }
        }

        if(maxRuns < V[0])
        {
            ++nu[0];
        }

        for(size_t j = 0; j <= K; ++j)
        {
            if(maxRuns == V[j])
            {
                ++nu[j];
            }
        }

        if(maxRuns > V[K])
        {
            ++nu[K];
        }
    }

    double chi2 = 0.0;
    for (size_t i = 0; i <= K; ++i )
    {
        chi2 += std::pow(static_cast<double>(nu[i]) - static_cast<double>(numberOfBlock) * pi[i], 2) / (static_cast<double>(numberOfBlock) * pi[i]);
    }

    double P = boost::math::gamma_q((static_cast<double>(K)/2.0), chi2/2.0);

    return {P >= threshold, P};
}

} // namespace nistpp
