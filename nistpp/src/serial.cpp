#include "nistpp/bits_storage.h"
#include "nistpp/types.h"

#include <cstddef>
#include <nistpp/tests.h>

#include <cmath>

#include <boost/math/special_functions/gamma.hpp>
#include <vector>

namespace nistpp
{

constexpr double  threshold     = 0.01;

double psi2(const BitsStorage::bits_t bits, std::size_t m)
{
    if(m == 0)
    {
        return 0.0;
    }

    const std::size_t powLen = std::pow(2, m+1) - 1;
    std::vector<std::size_t> P(powLen, 0);

    const auto n = bits.size();
    for(std::size_t i = 0; i < n; ++i)
    {
        std::size_t k = 1;
        for(std::size_t j = 0; j < m; ++j)
        {
            if(bits[(i + j) % n] == 0)
            {
                k *= 2;
            }
            else
            {
                k = k * 2 + 1;
            }
        }

        ++P[k-1];
    }

    double sum = 0;
    for(std::size_t i = std::pow(2, m) - 1; i < powLen; ++i)
    {
        sum += std::pow(P[i], 2);
    }

    sum = (sum * std::pow(2, m) / static_cast<double>(n)) - static_cast<double>(n);

    return sum;
}

return_t SerialTest(const BitsStorage& data, std::size_t M)
{
    const auto& bits = data.GetBits();

    const auto psim0 = psi2(bits, M);
    const auto psim1 = psi2(bits, M - 1);
    const auto psim2 = psi2(bits, M - 2);

    const auto del1 = psim0 - psim1;
    const auto del2 = psim0 - 2.0 * psim1 + psim2;

    const auto powM     = std::pow(2, M - 1) / 2.0;
    const auto p_value1 = boost::math::gamma_q(powM, del1/2.0);
    const auto p_value2 = boost::math::gamma_q(powM, del2/2.0);

    const auto minP = std::fmin(p_value1, p_value2);

    return {minP >= threshold, minP};
}

} // namespace nistpp
