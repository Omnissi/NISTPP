#include <cstddef>
#include <nistpp/tests.h>

#include <cmath>

#include "nistpp/types.h"
#include "sprout/math/sqrt.hpp"
#include "sprout/valarray/exponential.hpp"

namespace nistpp
{

constexpr double  threshold     = 0.01;

return_t UniversalTest(const BitsStorage& data)
{
    constexpr std::array<double, 17> expected_value = {{0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
                                            8.1764248, 9.1723243, 10.170032, 11.168765,
                                            12.168070, 13.167693, 14.167488, 15.167379 }};
    constexpr std::array<double, 17>  variance = {{ 0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
                                       3.401, 3.410, 3.416, 3.419, 3.421 }};

    std::size_t L = 5;
    const auto n = data.NumberOfBits();
    if ( n >= 387840 )     L = 6;
    if ( n >= 904960 )     L = 7;
    if ( n >= 2068480 )    L = 8;
    if ( n >= 4654080 )    L = 9;
    if ( n >= 10342400 )   L = 10;
    if ( n >= 22753280 )   L = 11;
    if ( n >= 49643520 )   L = 12;
    if ( n >= 107560960 )  L = 13;
    if ( n >= 231669760 )  L = 14;
    if ( n >= 496435200 )  L = 15;
    if ( n >= 1059061760 ) L = 16;

    const std::size_t Q = 10 * static_cast<std::size_t>(std::pow(2, L));
    const std::size_t K = static_cast<std::size_t>(std::floor(n/L)) - Q;
    const std::size_t p = static_cast<std::size_t>(std::pow(2, L));

    std::vector<std::size_t> T(p);

    const double c = 0.7 - 0.8/static_cast<double>(L)
            + (4.0 + 32.0/static_cast<double>(L)) * std::pow(K, -3.0/static_cast<double>(L))/15.0;
    const double sigma = c * std::sqrt(variance[L]/static_cast<double>(K));

    const auto& bits = data.GetBits();
    for(std::size_t i = 1; i <= Q; ++i)
    {
        std::size_t decRep = 0;
        for(std::size_t j = 0; j < L; ++j)
        {
            if(bits[(i - 1) * L + j])
            {
                decRep += static_cast<std::size_t>(std::pow(2, L - 1 - j));
            }
        }

        T[decRep] = i;
    }

    double sum = 0;
    for(std::size_t i = Q + 1; i <= Q + K; ++i)
    {
        std::size_t decRep = 0;
        for(std::size_t j = 0; j < L; ++j)
        {
            if(bits[(i - 1) * L + j])
            {
                decRep += static_cast<std::size_t>(std::pow(2, L - 1 - j));
            }
        }
        sum += std::log(i - T[decRep]) / sprout::log(2);
        T[decRep] = i;
    }

    const double phi     = sum / static_cast<double>(K);
    const double arg     = std::fabs(phi-expected_value[L])/(sprout::sqrt(2) * sigma);
    const auto   p_value = std::erfc(arg);

    return {p_value >= threshold, p_value};
}

} // namespace nistpp
