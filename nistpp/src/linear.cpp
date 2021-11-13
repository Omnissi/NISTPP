#include "nistpp/bits_storage.h"

#include <cstddef>
#include <nistpp/tests.h>

#include <cmath>

#include <boost/math/special_functions/gamma.hpp>
#include <vector>

namespace nistpp
{

return_t LinearComplexityTest(const BitsStorage& data, std::size_t M)
{
    constexpr std::size_t K     = 6;
    constexpr double      pi[7] = {0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833};

    const auto        N     = static_cast<size_t>(std::floor(data.NumberOfBits() / M));
    const int8_t      sign  = ((M + 1) % 2) == 0 ? -1 : 1;
    const double      mean  = static_cast<double >(M)/2.0 + (9.0 + sign) / 36.0 - 1.0 / std::pow(2, M) * (static_cast<double >(M) / 3.0 + 2.0/9.0);

    std::array<std::atomic_uint32_t, 7> nu = {};

    const auto& bits = data.GetBits();

#pragma omp parallel for
    for(std::size_t i = 0; i < N; ++i)
    {
        std::vector<BitsStorage::bits_t::value_type> T(M);
        std::vector<BitsStorage::bits_t::value_type> C(M);
        std::vector<BitsStorage::bits_t::value_type> B(M);

        C[0] = true;
        B[0] = true;

        std::size_t d;
        ssize_t     L = 0;
        ssize_t     m = -1;

        std::size_t N_ = 0;
        while(N_ < M)
        {
            std::size_t constInd = i * M + N_;
            d = static_cast<decltype(d)>(bits[constInd]);
            for(std::size_t j = 1; j <= static_cast<std::size_t>(L); ++j)
            {
                if(bits[constInd - j] && C[j])
                {
                    ++d;
                }
            }

            d %= 2;
            if(d)
            {
                std::copy(C.begin(), C.end(), T.begin());

                constInd = N_ - static_cast<std::size_t>(m);

                for(std::size_t j = 0; (j + constInd) < M; ++j)
                {
                    C[j + constInd] = (C[j + constInd] + B[j]) % 2;
                }

                if(L <= static_cast<decltype(L)>(N_ / 2))
                {
                    L = static_cast<decltype(L)>(N_) + 1 - L;
                    m = static_cast<decltype(m)>(N_);
                    std::copy(T.begin(), T.end(), B.begin());
                }
            }
            ++N_;
        }

        const double T_ = sign * (static_cast<double >(L) - mean) + 2.0/9.0;

        if (T_ <= -2.5)
        {
            ++nu[0];
        }
        else if (T_ > -2.5 && T_ <= -1.5)
        {
            ++nu[1];
        }
        else if (T_ > -1.5 && T_ <= -0.5)
        {
            ++nu[2];
        }
        else if (T_ > -0.5 && T_ <= 0.5)
        {
            ++nu[3];
        }
        else if (T_ >  0.5 && T_ <= 1.5)
        {
            ++nu[4];
        }
        else if (T_ >  1.5 && T_ <= 2.5)
        {
            ++nu[5];
        }
        else
        {
            ++nu[6];
        }
    }

    double chi2 = 0;
    for(std::size_t i = 0; i < K + 1; ++i)
    {
        chi2 += std::pow(nu[i] - static_cast<double >(N) * pi[i], 2) / (static_cast<double >(N) * pi[i]);
    }

    const double p_value = boost::math::gamma_q(K/2.0, chi2/2.0);

    return {p_value >= threshold, p_value};
}

} // namespace nistpp
