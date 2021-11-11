#include "nistpp/bits_storage.h"
#include "nistpp/types.h"

#include <nistpp/tests.h>

#include <cmath>
#include <vector>

#include <boost/math/distributions.hpp>

#include <sprout/math.hpp>

namespace nistpp
{

constexpr double  threshold     = 0.01;

double normal(double x)
{
    constexpr auto sqrt2 = sprout::sqrt(2);

    return (1.0 + std::erf(x/sqrt2)) / 2.0;
}

return_t CumulativeSumsTest(const BitsStorage& data)
{
    int32_t S   = 0;
    int32_t sup = 0;
    int32_t inf = 0;
    int32_t z;
    int32_t zrev;

    const auto& bits = data.GetBits();

    // @todo Проверить, что за херня тут происходит
    for(const auto& el : bits)
    {
        el ? ++S : --S;
        if(S > sup)
        {
            ++sup;
        }
        if(S < inf)
        {
            --inf;
        }

        z = (sup > -inf) ? sup : -inf;
        zrev = (sup - S > S - inf) ? sup - S : S - inf;
    }

    const int32_t n     = static_cast<int32_t>(bits.size());
    const double  sqrtn = std::sqrt(n);

    int32_t begin = (-n / z + 1) / 4;
    int32_t end   = (n / z - 1) / 4;

    double sum1 = 0;
    for(auto k = begin; k <= end; ++k)
    {
        sum1 += normal(((4 * k + 1) * z) / sqrtn);
        sum1 -= normal(((4 * k - 1) * z) / sqrtn);
    }

    double sum2 = 0;
    for(auto k = (-n / z - 3) / 4; k <= end; ++k)
    {
        sum2 += normal(((4 * k + 3) * z) / sqrtn);
        sum2 -= normal(((4 * k + 1) * z) / sqrtn);
    }

    const double p_value_f = 1.0 - sum1 + sum2;

    begin = (-n / zrev + 1) / 4;
    end   = (n / zrev - 1) / 4;

    sum1 = 0;
    for(auto k = begin; k <= end; ++k)
    {
        sum1 += normal(((4 * k + 1) * zrev) / sqrtn);
        sum1 -= normal(((4 * k - 1) * zrev) / sqrtn);
    }

    sum2 = 0;
    for(auto k = (-n / zrev - 3) / 4; k <= end; ++k)
    {
        sum2 += normal(((4 * k + 3) * zrev) / sqrtn);
        sum2 -= normal(((4 * k + 1) * zrev) / sqrtn);
    }

    const double p_value_b = 1.0 - sum1 + sum2;
    const double minP      = std::fmin(p_value_b, p_value_f);

    return {minP >= threshold, minP};
}

} // namespace nistpp