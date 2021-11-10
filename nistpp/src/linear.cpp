#include "nistpp/bits_storage.h"

#include <cstddef>
#include <nistpp/tests.h>

#include <cmath>

#include <boost/math/special_functions/gamma.hpp>
#include <vector>

namespace nistpp
{

constexpr double  threshold     = 0.01;

return_t LinearComplexityTest(const BitsStorage& data, std::size_t M)
{
    constexpr std::size_t K     = 6;
    constexpr double      pi[7] = {0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833};

    const std::size_t N     = std::floor(data.NumberOfBits() / M);
    const int8_t      sign  = ((M + 1) % 2) == 0 ? -1 : 1;
    const double      mean  = M/2.0 + (9.0 + sign) / 36.0 - 1.0 / std::pow(2, M) * (M / 3.0 + 2.0/9.0);

    std::vector<BitsStorage::bits_t::value_type> T(M);
    std::vector<BitsStorage::bits_t::value_type> C(M);
    std::vector<BitsStorage::bits_t::value_type> P(M);
    std::vector<BitsStorage::bits_t::value_type> B(M);

    std::array<double, 7> nu = {};

    const auto& bits = data.GetBits();

#define fill_zero(arr) std::fill(arr.begin(), arr.end(), 0)
    for(std::size_t i = 0; i < N; ++i)
    {
        fill_zero(T);
        fill_zero(C);
        fill_zero(P);
        fill_zero(B);

        C[0] = 1;
        B[0] = 1;

        int32_t     L = 0;
        std::size_t d = 0;
        int32_t     m = -1;

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
                // This solution is slow
                //std::copy(C.begin(), C.end(), T.begin());
                //fill_zero(P);
                // This is faster
                for(std::size_t j = 0; j < M; ++j)
                {
                    T[j] = C[j];
                    P[j] = 0;
                }

                constInd = N_ - m;
                for(std::size_t j = 0; j < M; ++j)
                {
                    if(B[j] == 1)
                    {
                        P[j + constInd] = 1;
                    }
                }

                for(std::size_t j = 0; j < M; ++j)
                {
                    C[j] = (C[j] + P[j]) % 2;
                }

                if(L <= static_cast<decltype(L)>(N_ / 2))
                {
                    L = N_ + 1 - L;
                    m = N_;
                    std::copy(T.begin(), T.end(), B.begin());
                }
            }
            ++N_;
        }

        const double T_ = sign * (L - mean) + 2.0/9.0;

        if (T_ <= -2.5)
            ++nu[0];
        else if (T_ > -2.5 && T_ <= -1.5)
            ++nu[1];
        else if (T_ > -1.5 && T_ <= -0.5)
            ++nu[2];
        else if (T_ > -0.5 && T_ <= 0.5)
            ++nu[3];
        else if (T_ >  0.5 && T_ <= 1.5)
            ++nu[4];
        else if (T_ >  1.5 && T_ <= 2.5)
            ++nu[5];
        else
            ++nu[6];
    }
#undef fill_zero

    double chi2 = 0;
    for(std::size_t i = 0; i < K + 1; ++i)
    {
        chi2 += std::pow(nu[i] - N * pi[i], 2) / (N * pi[i]);
    }

    const double p_value = boost::math::gamma_q(K/2.0, chi2/2.0);

    return {p_value >= threshold, p_value};
}

} // namespace nistpp
