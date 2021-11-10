#include <algorithm>
#include <nistpp/tests.h>
#include <nistpp/math_helpers.h>

#include <sprout/math.hpp>

#include "templates/templates.h"

#include "templates/template2.hpp"
#include "templates/template3.hpp"
#include "templates/template4.hpp"
#include "templates/template5.hpp"
#include "templates/template6.hpp"
#include "templates/template7.hpp"
#include "templates/template8.hpp"
#include "templates/template9.hpp"
#include "templates/template10.hpp"
#include "templates/template11.hpp"
#include "templates/template12.hpp"
#include "templates/template13.hpp"
#include "templates/template14.hpp"
#include "templates/template15.hpp"
//#include "templates/template16.hpp"

#include <cmath>

#include <boost/math/special_functions/gamma.hpp>
#include <utility>

namespace nistpp
{

constexpr double  threshold             = 0.01;
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


    const double lambda   = (M-m+1)/std::pow(2, m);
    const double varWj    = M*(1.0/std::pow(2.0, m) - (2.0*m-1.0)/std::pow(2.0, 2.0*m));
    const double sqrVarWj = std::pow(varWj, 0.5);

    auto        numberOfRows = GetNumberOfRows(m);

    double minP         = std::numeric_limits<double>::max();
    const auto& bits    = data.GetBits();
    numberOfRows        = std::min(numberOfRows, maxNumOfTemplates);
    P.resize(numberOfRows);

#pragma omp parallel for
    for(std::size_t i = 0; i < numberOfRows; ++i)
    {
        std::vector<uint32_t> Wj(N);
        auto begin  = bits.begin();
        auto end    = begin + m;

        auto sec_it = GetTemplatesSequence(i, m);
        auto tbegin = std::get<0>(sec_it);
        auto tend   = std::get<1>(sec_it);

        for(std::size_t j = 0; j < N; ++j)
        {
            std::size_t W_obs = 0;
            for(std::size_t k = 0; k < M-m+1; ++k, ++begin, ++end)
            {
                if(std::equal(tbegin, tend, begin, end))
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
            chi2 += std::pow((static_cast<double>(Wj[j]) - lambda)/sqrVarWj, 2);
        }

        P[i] = boost::math::gamma_q(N/2.0, chi2/2.0);
        minP = std::fmin(P[i], minP);
    }

    return {minP >= threshold, minP};
}

} // namespace nistpp
