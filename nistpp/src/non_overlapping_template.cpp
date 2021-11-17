#include <algorithm>
#include <nistpp/tests.h>
#include <nistpp/math_helpers.h>

#include "templates/templates.h"

#include <cmath>

#include <utility>

namespace nistpp
{

constexpr std::size_t maxNumOfTemplates = 148;

return_t NonOverlappingTemplateTest(const BitsStorage& data, std::size_t m, std::vector<double> &P)
{
    if( 2 > m || m > 16)
    {
        throw std::logic_error("In this implementation 2 <= m <= 16 !");
    }

    constexpr size_t N = 8;
    const auto numberOfBits = data.NumberOfBits();
    const size_t M = numberOfBits/N;


    const double lambda   = static_cast<double>(M-m+1)/std::pow(2, m);
    const double varWj    = static_cast<double>(M)*(1.0/std::pow(2.0, m) - (2.0*static_cast<double>(m)-1.0)/std::pow(2.0, 2.0*static_cast<double >(m)));
    const double sqrVarWj = std::pow(varWj, 0.5);

    auto        numberOfRows = GetNumberOfRows(m);

    const auto& bits    = data.GetBits();
    numberOfRows        = std::min(numberOfRows, maxNumOfTemplates);
    P.resize(numberOfRows);

#pragma omp parallel for default(shared)
    for(std::size_t i = 0; i < numberOfRows; ++i)
    {
        std::vector<std::size_t> Wj(N);

        auto sec_it = GetTemplatesSequence(i, m);
        auto tbegin = std::get<0>(sec_it);
        auto tend   = std::get<1>(sec_it);

        for(std::size_t j = 0; j < N; ++j)
        {
            auto begin  = bits.begin() + static_cast<ssize_t>(j * M);
            auto end    = begin + static_cast<ssize_t>(m);

            std::size_t W_obs = 0;
            for(std::size_t k = 0; k < M-m+1; ++k, ++begin, ++end)
            {
                if(std::equal(tbegin, tend, begin, end))
                {
                    ++W_obs;
                    begin += static_cast<ssize_t>(m - 1);
                    end += static_cast<ssize_t>(m - 1);
                    k += m-1;
                }
            }
            Wj[j] = W_obs;
        }

        double chi2 = 0;
        for (std::size_t j = 0; j < N; ++j )
        {
            chi2 += std::pow((static_cast<double>(Wj[j]) - lambda)/sqrVarWj, 2);
        }

        P[i] = igamc(N/2.0, chi2/2.0);
    }

    const double minP = *std::min_element(P.begin(), P.end());

    return {minP >= threshold, minP};
}

} // namespace nistpp
