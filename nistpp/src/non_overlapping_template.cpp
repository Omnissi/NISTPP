#include <nistpp/tests.h>
#include <nistpp/math_helpers.h>

#include <sprout/math.hpp>

#include "templates/templates.h"

#include <cmath>

#include <iostream>

#include <boost/math/special_functions/gamma.hpp>

namespace nistpp
{

constexpr double  threshold     = 0.01;
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

    std::vector<uint32_t> Wj(N);

    const double lambda   = (M-m+1)/pow(2, m);
    const double varWj    = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
    const double sqrVarWj = std::pow(varWj, 0.5);

    std::size_t numberOfRows = 0;
    auto funcTempl = GetTemplatesFunction(m, numberOfRows);

    P.resize(numberOfRows);
    double minP = std::numeric_limits<double>::max();
    auto bits   = data.GetBits();
    for(std::size_t i = 0; i < std::min(numberOfRows, maxNumOfTemplates) ; ++i)
    {
        auto begin  = bits.begin();
        auto end    = begin + m;

        for(std::size_t j = 0; j < N; ++j)
        {
            std::size_t W_obs = 0;
            for(std::size_t k = 0; k < M-m+1; ++k, ++begin, ++end)
            {
                if(funcTempl(i, begin, end))
                {
                    ++W_obs;
                    begin += m-1;
                    end += m-1;
                    k += m-1;
                }
            }
            Wj[j] = W_obs;
        }

        double chi2 = 0;
        for (std::size_t j = 0; j < N; ++j )
        {
            chi2 += pow((static_cast<double>(Wj[j]) - lambda)/sqrVarWj, 2);
        }

        P[i] = boost::math::gamma_q(N/2.0, chi2/2.0);
        minP = std::fmin(P[i], minP);
    }

    return {minP >= threshold, minP};
}

} // namespace nistpp
