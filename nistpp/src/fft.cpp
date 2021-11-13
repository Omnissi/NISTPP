#include <nistpp/tests.h>

#include <kissfft/kissfft.hh>

#include <valarray>
#include <cmath>

namespace nistpp
{

return_t FftTest(const BitsStorage &data)
{
    const auto N = data.NumberOfBits();

    using FFT   = kissfft<double>;
    using cpx_t = std::complex<double>;

    std::vector<cpx_t> x(N);
    for(std::size_t i = 0; i < N; ++i)
    {
        x[i] = data[i] ? 1 : -1;
    }

    std::vector<cpx_t> res(N);

    FFT fft(N, false);
    fft.transform(x.data(), res.data());

    std::vector<double> m(N/2 + 1);
    for(std::size_t i = 0; i < m.size(); ++i)
    {
        m[i] = static_cast<double>(std::abs(res[i]));
    }

    size_t count = 0;
    const auto upperBound = std::sqrt(2.995732274 * static_cast<double>(N));
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
