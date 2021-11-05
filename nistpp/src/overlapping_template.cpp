#include <cstddef>
#include <nistpp/tests.h>
#include <nistpp/math_helpers.h>

#include <sprout/math.hpp>

#include "templates/templates.h"

#include <cmath>

#include <boost/math/special_functions/gamma.hpp>

namespace nistpp
{

constexpr double  threshold             = 0.01;

double computePi(std::size_t u, double eta)
{
    if (u == 0)
    {
        return std::exp(-eta);
    }

    double sum = 0.0;
    for(std::size_t l = 1; l<=u; ++l)
    {
        sum += exp(-eta - u * std::log(2) + l * std::log(eta) - std::lgamma(l + 1)
                   + std::lgamma(u) - std::lgamma(l) - std::lgamma(u - l + 1));
    }

    return sum;
}

return_t OverlappingTemplateTest(const BitsStorage& data, std::size_t m)
{
    constexpr std::size_t K = 5;
    constexpr std::size_t M = 1032;
    const     std::size_t N = data.NumberOfBits() / M;

    std::array<double, 6> pi = {};
    {
        double lambda = (M - m + 1) / std::pow(2, m);
        double eta    = lambda / 2.0;
        double sum    = 0.0;

        for(std::size_t i = 0; i < K; ++i)
        {
            pi[i] = computePi(i, eta);
            sum += pi[i];
        }

        pi[K] = 1 - sum;
    }

    const std::vector<bool> test_seq(m, true);
    const auto& eps     = data.GetBits();
    auto        it_eps  = eps.begin();

    std::array<std::size_t, 6> nu = {};
    for(std::size_t i = 0; i < N; ++i)
    {
        std::size_t W_obs = 0;
        for(std::size_t j = 0; j < M - m + 1; ++j, ++it_eps)
        {
            if(std::equal(test_seq.begin(), test_seq.end(), it_eps))
            {
                ++W_obs;
            }
        }

        it_eps += m - 1;
        W_obs   = std::min(W_obs, nu.size() - 1);
        ++nu[W_obs];
    }

    double chi2 = 0.0;
    for(std::size_t i = 0; i < K+1; ++i)
    {
        chi2 += pow(static_cast<double>(nu[i]) - static_cast<double>(N*pi[i]), 2)/(static_cast<double>(N*pi[i]));
    }

    double p_value = boost::math::gamma_q(K/2.0, chi2/2.0);

    return {p_value >= threshold, p_value};
}

} // namespace nistpp
