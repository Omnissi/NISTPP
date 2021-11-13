#include <nistpp/tests.h>

#include "dfft.h"

#include <valarray>
#include <cmath>

namespace nistpp
{

return_t FftTest(const BitsStorage &data)
{
    const auto N = data.NumberOfBits();

//    CArray x(N);
    std::vector<double> x(N);
    for(std::size_t i = 0; i < N; ++i)
    {
        x[i] = data[i] ? 1 : -1;
    }

    std::vector<double> wsave(2 * N);
    int ifac[15] = {};
    __ogg_fdrffti(static_cast<int>(N), wsave.data(), ifac);		/* INITIALIZE WORK ARRAYS */
    __ogg_fdrfftf(static_cast<int>(N), x.data(), wsave.data(), ifac);	/* APPLY FORWARD FFT */
//    fft(x);

    std::valarray<double> m(N/2 + 1);
    m[0] = std::fabs(x[0]);

    for(size_t i = 0; i < N/2; ++i)
    {
        m[i+1] = /*std::fabs(x[i+1]);*/sqrt(pow(x[2*i+1],2)+pow(x[2*i+2],2));
    }
    size_t count = 0;
    auto upperBound = std::sqrt(2.995732274 * static_cast<double>(N));
    for (size_t i = 0; i < N/2; ++i)
    {
        if (m[i] < upperBound)
        {
            ++count;
        }
    }

    double N_l = static_cast<double>(count);       /* number of peaks less than h = sqrt(3*n) */
    double N_o = static_cast<double>(0.95*static_cast<double>(N)/2.0);
    double d   = (N_l - N_o)/std::sqrt(static_cast<double>(N)/4.0*0.95*0.05);
    double P   = std::erfc(std::fabs(d)/std::sqrt(2.0));

    return {P >= threshold, P};
}

} // namespace nistpp
