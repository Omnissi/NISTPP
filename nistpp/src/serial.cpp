#include "nistpp/bits_storage.h"
#include "nistpp/types.h"
#include <nistpp/tests.h>
#include <nistpp/math_helpers.h>

#include <cstddef>
#include <cmath>
#include <vector>

namespace nistpp
{

double psi2(const BitsStorage::bits_t& bits, std::size_t m)
{
    if(m == 0)
    {
        return 0.0;
    }

    const std::size_t powLen = static_cast<std::size_t>(std::pow(2, m+1)) - 1;
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
    for(std::size_t i = static_cast<std::size_t>(std::pow(2, m)) - 1; i < powLen; ++i)
    {
        sum += std::pow(P[i], 2);
    }

    sum = (sum * std::pow(2, m) / static_cast<double>(n)) - static_cast<double>(n);

    return sum;
}

return_t SerialTest(const BitsStorage &data, std::size_t M, std::array<double, 2> &P)
{
    const auto& bits = data.GetBits();

    const auto psim0 = psi2(bits, M);
    const auto psim1 = psi2(bits, M - 1);
    const auto psim2 = psi2(bits, M - 2);

    const auto del1 = psim0 - psim1;
    const auto del2 = psim0 - 2.0 * psim1 + psim2;

    P[0] = igamc(std::pow(2, M - 1) / 2.0, del1/2.0);
    P[1] = igamc(std::pow(2, M - 2) / 2.0, del2/2.0);

    const auto minP = std::fmin(P[0], P[1]);

    return {minP >= threshold, minP};
}

} // namespace nistpp
