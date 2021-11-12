#include "nistpp/bits_storage.h"
#include "nistpp/types.h"

#include <nistpp/tests.h>

#include <cmath>
#include <vector>

#include <boost/math/special_functions/gamma.hpp>

namespace nistpp
{

constexpr double threshold = 0.01;

return_t RandomExcursionsTest(const BitsStorage &data, std::array<double, 8>& P)
{
    constexpr std::array<ssize_t, 8> stateX = {{-4, -3, -2, -1, 1, 2, 3, 4}};
    constexpr std::array<std::array<double, 6>, 5> pi = {{{0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000},
                                                                 {0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000},
                                                                 {0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625},
                                                                 {0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143},
                                                                 {0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051}}};

    const auto maxIterations = std::max(static_cast<std::size_t>(1000), data.NumberOfBits() / 100);
    std::vector<ssize_t>     S_k(data.NumberOfBits());
    std::vector<std::size_t> cycle(maxIterations);

    const auto& bits = data.GetBits();

    std::size_t J = 0;
    S_k[0] = 2 * bits[0] - 1;

    for(std::size_t i = 1; i < bits.size(); ++i)
    {
        S_k[i] = S_k[i - 1] + 2 * static_cast<ssize_t>(bits[i]) - 1;
        if(S_k[i] == 0)
        {
            ++J;
            if (J > maxIterations)
            {
                throw std::runtime_error(std::string("In ")
                                         + std::string(__func__) +
                                         std::string(":  EXCEEDING THE MAX NUMBER OF CYCLES EXPECTED"));
            }

            cycle[J] = i;
        }
    }

    if(S_k[bits.size() - 1] != 0)
    {
        ++J;
    }

    cycle[J] = bits.size();

    if(static_cast<double>(J) < std::fmax(0.005 * std::pow(bits.size(), 0.5), 500.0))
    {
        throw std::runtime_error("WARNING:  TEST NOT APPLICABLE.  THERE ARE AN INSUFFICIENT NUMBER OF CYCLES.");
    }

    std::array<std::array<std::size_t, 8>, 6> nu      = {};
    std::array<std::size_t, 8>                counter = {};

    std::size_t cycleStart = 0;
    std::size_t cycleStop  = cycle[1];

    for(std::size_t j = 1; j <= J; ++j)
    {
        std::fill(counter.begin(), counter.end(), 0);
        for(std::size_t i = cycleStart; i < cycleStop; ++i)
        {
            if((S_k[i] >= 1 && S_k[i] <= 4) || (S_k[i] >= -4 && S_k[i] <= -1))
            {
                if(S_k[i] < 0)
                {
                    ++counter[static_cast<std::size_t>(S_k[i] + 4)];
                }
                else
                {
                    ++counter[static_cast<std::size_t>(S_k[i] + 3)];
                }
            }
        }

        cycleStart = cycle[j]+1;
        if(j < J)
        {
            cycleStop = cycle[j+1];
        }

        for(std::size_t i = 0; i < counter.size(); ++i)
        {
            const auto index = std::min(counter[i], static_cast<std::size_t>(5));
            ++nu[index][i];
        }
    }

    double minP = std::numeric_limits<double>::max();
    for(std::size_t i = 0; i < stateX.size(); ++i)
    {
        const auto x   = stateX[i];
        double     sum = 0;
        for(std::size_t j = 0; j < nu.size(); ++j)
        {
            const auto tmp = static_cast<double>(J) * pi[static_cast<std::size_t>(std::abs(x))][j];
            sum += std::pow(static_cast<double>(nu[j][i]) - tmp, 2) / tmp;
        }

        P[i] = boost::math::gamma_q(2.5, sum / 2.0);
        minP = std::fmin(minP, P[i]);
    }

    return {minP >= threshold, minP};
}

} // namespace nistpp
