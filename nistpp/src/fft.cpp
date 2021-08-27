#include <nistpp/tests.h>

#include <nistpp/dfft.h>

#include <valarray>
#include <complex>
#include <cmath>

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

void fft(CArray &x)
{
    const size_t N = x.size();
    if (N <= 1) return;

    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];

    // conquer
    fft(even);
    fft(odd);

    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

return_t FftTest(const BitsStorage &data)
{
    const auto N = data.NumberOfBits();

//    CArray x(N);
    std::vector<double> x(N);
    const auto& bits = data.GetBits();
    for(std::size_t i = 0; i < bits.size(); ++i)
    {
        const auto& word = bits[i];
        const auto ind_i = i * BitsStorage::numberOfBitsInWord;
        for(std::size_t j = 0; j < BitsStorage::numberOfBitsInWord; ++j)
        {
            x[ind_i + j] = word[BitsStorage::numberOfBitsInWord - 1 - j] ? 1 : -1;//Complex(1, 0) : Complex(0, -1);
        }
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
        m[i+1] = sqrt(pow(x[2*i+1],2)+pow(x[2*i+2],2));
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
    double percentile = static_cast<double>(count)/(N/2.0)*100.0;
    double N_l = static_cast<double>(count);       /* number of peaks less than h = sqrt(3*n) */
    double N_o = static_cast<double>(0.95*N/2.0);
    double d   = (N_l - N_o)/std::sqrt(N/4.0*0.95*0.05);
    double P   = std::erfc(fabs(d)/std::sqrt(2.0));

    return {P >= threshold, P};
}

} // namespace nistpp
