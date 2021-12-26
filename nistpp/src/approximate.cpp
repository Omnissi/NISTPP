#include "nistpp/bits_storage.h"
#include "nistpp/types.h"
#include "nistpp/math_helpers.h"

#include <nistpp/tests.h>

#include <cmath>
#include <vector>

#include <gcem.hpp>

namespace nistpp
{

return_t ApproximateEntropyTest(const BitsStorage& data, std::size_t M)
{
    std::array<double, 2> ApEn = {};
    std::size_t r = 0;

    const auto& bits = data.GetBits();

    for(std::size_t blockSize = M; blockSize <= M + 1; ++blockSize)
    {
        if(blockSize == 0)
        {
            ApEn[0] = 0.0;
            ++r;
            continue;
        }

        const std::size_t powLen = static_cast<std::size_t>(std::pow(2, blockSize + 1)) - 1;
        std::vector<std::size_t> P(powLen);

        for(std::size_t i = 0; i < bits.size(); ++i)
        {
            std::size_t k = 1;
            for(std::size_t j = 0; j < blockSize; ++j)
            {
                k <<= 1;
                if(bits[(i + j) % bits.size()])
                {
                    ++k;
                }
            }
            ++P[k - 1];
        }

        double sum = 0;
        const std::size_t powConst = static_cast<std::size_t>(std::pow(2, blockSize));
        std::size_t index =  powConst - 1;
        for(std::size_t i = 0; i < powConst; ++i, ++index)
        {
            if(P[index] > 0)
            {
                sum += static_cast<double>(P[index]) * std::log(static_cast<double>(P[index]) / static_cast<double>(bits.size()));
            }
        }

        sum /= static_cast<double>(bits.size());
        ApEn[r] = sum;
        ++r;
    }

    const auto apen    = ApEn[0] - ApEn[1];
    const auto chi2    = 2.0 * static_cast<double>(bits.size()) * (gcem::log(2) - apen);
    const auto p_value = igamc(std::pow(2, M-1), chi2/2.0);

    return {p_value >= threshold, p_value};
}

} // namespace nistpp
