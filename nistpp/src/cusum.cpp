#include "nistpp/bits_storage.h"
#include "nistpp/types.h"

#include <nistpp/tests.h>

#include <cmath>

#include <boost/math/distributions.hpp>

#include <sprout/math.hpp>

namespace nistpp
{

double normal(double x)
{
    constexpr auto sqrt2 = sprout::sqrt(2);

    return (1.0 + std::erf(x/sqrt2)) / 2.0;
}

return_t CumulativeSumsTest(const BitsStorage &data, std::array<double, 2> &P)
{
    ssize_t S    = 0;
    ssize_t sup  = 0;
    ssize_t inf  = 0;
    ssize_t z    = 0;
    ssize_t zrev = 0;

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

    const ssize_t n     = static_cast<ssize_t>(bits.size());
    const double  sqrtn = std::sqrt(n);

    ssize_t begin = (-n / z + 1) / 4;
    ssize_t end   = (n / z - 1) / 4;

    double sum1 = 0;
    for(auto k = begin; k <= end; ++k)
    {
        sum1 += normal((static_cast<double>((4 * k + 1) * z)) / sqrtn);
        sum1 -= normal((static_cast<double>((4 * k - 1) * z)) / sqrtn);
    }

    double sum2 = 0;
    for(auto k = (-n / z - 3) / 4; k <= end; ++k)
    {
        sum2 += normal((static_cast<double>((4 * k + 3) * z)) / sqrtn);
        sum2 -= normal((static_cast<double>((4 * k + 1) * z)) / sqrtn);
    }

    P[0] = 1.0 - sum1 + sum2;

    begin = (-n / zrev + 1) / 4;
    end   = (n / zrev - 1) / 4;

    sum1 = 0;
    for(auto k = begin; k <= end; ++k)
    {
        sum1 += normal((static_cast<double>((4 * k + 1) * zrev)) / sqrtn);
        sum1 -= normal((static_cast<double>((4 * k - 1) * zrev)) / sqrtn);
    }

    sum2 = 0;
    for(auto k = (-n / zrev - 3) / 4; k <= end; ++k)
    {
        sum2 += normal((static_cast<double>((4 * k + 3) * zrev)) / sqrtn);
        sum2 -= normal((static_cast<double>((4 * k + 1) * zrev)) / sqrtn);
    }

    P[1] = 1.0 - sum1 + sum2;
    const double minP      = std::fmin(P[0], P[1]);

    return {minP >= threshold, minP};
}

} // namespace nistpp