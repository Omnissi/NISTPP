#include <nistpp/tests.h>

#include <cmath>

#include <gcem.hpp>

namespace nistpp
{

return_t FrequencyTest(const BitsStorage &data)
{
    const auto numberOfBits = data.NumberOfBits();
    const auto ones         = data.NumberOfOnes();

    ssize_t Sn = static_cast<ssize_t>(2 * ones - (numberOfBits));
    double Sobs = static_cast<double>(std::abs(Sn)) / std::sqrt(numberOfBits);
    double P    = std::erfc(Sobs/gcem::sqrt(2.0));

    return {P >= threshold, P};
}

} // namespace nistpp
