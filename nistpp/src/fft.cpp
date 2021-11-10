#include <nistpp/tests.h>

#include "dfft.h"

#include <valarray>
#include <complex>
#include <cmath>
#include <unordered_map>

namespace nistpp
{

constexpr double  threshold     = 0.01;

using Complex   = std::complex<double>;
using CArray    = std::valarray<Complex>;

int reverse(int num, int lg_n) {
    int res = 0;
    for (int i = 0; i < lg_n; i++) {
        if (num & (1 << i))
            res |= 1 << (lg_n - 1 - i);
    }
    return res;
}

CArray fft(CArray &x)
{
    std::unordered_map<size_t, Complex> cash;

    size_t counter = 0;
    const double tmp = 2 * M_PI / static_cast<double>(x.size());
    CArray res(x.size());
    auto calcValue = [&](size_t k, size_t n) -> Complex
    {
        const size_t mp = k * n;
        auto it = cash.find(mp);
        if(it == cash.end())
        {
            auto res = Complex(std::cos(tmp * mp), std::sin(tmp * mp));
            cash[mp] = res;
            return res;
        }
        else
        {
            ++counter;
            return it->second;
        }
    };

    const auto N1 = x.size() - 1;
    for(size_t k = 0; k < N1; ++k)
    {
        for(size_t n = 0; n < N1; ++n)
        {
            res[k] += x[n] * calcValue(k, n);
        }
    }

    return res;
}

return_t FftTest(const BitsStorage &data)
{
    const auto N = data.NumberOfBits();

//    CArray x(N);
    std::vector<double> x(N);
    for(std::size_t i = 0; i < N; ++i)
    {
        x[i] = data[i] ? 1 : -1;
    }

    std::vector<double> wsave(2*N);
    int ifac[15] = {};
    __ogg_fdrffti(N, wsave.data(), ifac);		/* INITIALIZE WORK ARRAYS */
    __ogg_fdrfftf(N, x.data(), wsave.data(), ifac);	/* APPLY FORWARD FFT */
//    fft(x);

    std::valarray<double> m(N/2 + 1);
    m[0] = std::fabs(x[0]);

    for(size_t i = 0; i < N/2; ++i)
    {
        m[i+1] = /*std::fabs(x[i+1]);*/sqrt(pow(x[2*i+1],2)+pow(x[2*i+2],2));
    }
    size_t count = 0;
    auto upperBound = std::sqrt(2.995732274*N);
    for (size_t i = 0; i < N/2; ++i)
    {
        if (m[i] < upperBound)
        {
            ++count;
        }
    }

    double N_l = static_cast<double>(count);       /* number of peaks less than h = sqrt(3*n) */
    double N_o = static_cast<double>(0.95*N/2.0);
    double d   = (N_l - N_o)/std::sqrt(N/4.0*0.95*0.05);
    double P   = std::erfc(fabs(d)/std::sqrt(2.0));

    return {P >= threshold, P};
}

} // namespace nistpp
