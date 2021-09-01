#include <nistpp/tests.h>
#include <nistpp/math_helpers.h>

#include <cmath>

#include <sprout/math.hpp>

namespace nistpp
{

constexpr double  threshold     = 0.01;

return_t FrequencyTest(const BitsStorage &data)
{
    const auto& bits        = data.GetBits();
    const auto numberOfBits = data.NumberOfBits();
    const auto ones         = data.NumberOfOnes();

    int32_t Sn = static_cast<int32_t>(2 * ones - (numberOfBits));
    double Sobs = static_cast<double>(std::abs(Sn)) / std::sqrt(numberOfBits);
    double P    = std::erfc(Sobs/sprout::sqrt(2.0));

    return {P >= threshold, P};
}

} // namespace nistpp
