#include "nistpp/bits_storage.h"
#include "nistpp/types.h"

#include <nistpp/tests.h>

#include <cmath>
#include <vector>

#include <boost/math/special_functions/gamma.hpp>

namespace nistpp
{

return_t RandomExcursionsVariantTest(const BitsStorage& data, std::array<double, 18>& P)
{
    constexpr std::array<ssize_t, 18> stateX = {{-9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9}};

    const auto& bits = data.GetBits();
    std::vector<ssize_t> S_k(bits.size());

    std::size_t J = 0;
    S_k[0] = 2 * bits[0] - 1;
    for(std::size_t i = 1; i < bits.size(); ++i)
    {
        S_k[i] = S_k[i - 1] + 2 * bits[i] - 1;
        if(S_k[i] == 0)
        {
            ++J;
        }
    }

    if(S_k[bits.size() - 1] != 0)
    {
        ++J;
    }

    if(static_cast<double>(J) < std::fmax(0.005 * std::pow(bits.size(), 0.5), 500.0))
    {
        throw std::runtime_error("WARNING:  TEST NOT APPLICABLE.  THERE ARE AN INSUFFICIENT NUMBER OF CYCLES.");
    }

    double minP = std::numeric_limits<double>::max();
    for(std::size_t i = 0; i < stateX.size(); ++i)
    {
        const auto x = stateX[i];
        const auto count = std::count(S_k.begin(), S_k.end(), x);

        const auto arg = static_cast<double>(std::fabs(count - static_cast<ssize_t>(J)))
                / (std::sqrt(2.0 * static_cast<double>(J) * (4.0 * static_cast<double>(std::abs(x)) - 2)));
        P[i] = std::erfc(arg);
        minP = std::fmin(P[i], minP);
    }

    return {minP >= threshold, minP};
}

} // namespace nistpp
