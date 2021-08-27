#include <nistpp/tests.h>

#include <cmath>
#include <numeric>
#include <stdexcept>

#include "help_function.h"

#include <boost/math/special_functions/gamma.hpp>

namespace nistpp
{

constexpr double  threshold     = 0.01;

return_t BlockFrequencyTest(const BitsStorage &data, std::size_t M)
{
    const auto numberOfBits = data.NumberOfBits();
    if(numberOfBits < M)
    {
        throw std::invalid_argument("Size of block must be lower number of bits");
    }

    size_t numberOfBlocks = numberOfBits / M;

    double sum = 0.0;
    for(size_t i = 0; i < numberOfBlocks; ++i)
    {
        size_t blockSum = 0;
        for(size_t j = 0; j < M; ++j)
        {
            blockSum += data[i*M + j];
        }
        double v = (static_cast<double>(blockSum) / static_cast<double>(M)) - 0.5;
        sum += std::pow(v, 2);
    }

    double chi_squared = 4.0 * M * sum;
    double P           = boost::math::gamma_q(numberOfBlocks/2.0, chi_squared/2.0);

    return {P >= threshold, P};
}

} // namespace nistpp
