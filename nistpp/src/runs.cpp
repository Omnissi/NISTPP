#include <nistpp/tests.h>

#include <cmath>

#include "help_function.h"

namespace nistpp
{

return_t RunsTest(const BitsStorage& data)
{
    const auto numberOfBits = data.NumberOfBits();

    double pi = static_cast<double>(data.NumberOfOnes()) / static_cast<double>(numberOfBits);
    if(std::fabs(pi - 0.5) > (2.0 / std::sqrt(numberOfBits)))
    {
        return {false, 0.0};
    }

    size_t V = 1;
    for(size_t i = 1; i < numberOfBits; ++i)
    {
        if(data[i] != data[i-1])
        {
            ++V;
        }
    }

    double erfc_arg = fabs(static_cast<double>(V) - 2.0 * static_cast<double>(numberOfBits) * pi * (1-pi)) / (2.0 * pi * (1-pi) * std::sqrt(2*numberOfBits));
    double P        = std::erfc(erfc_arg);

    return {P >= threshold, P};
}

} // namespace nistpp
