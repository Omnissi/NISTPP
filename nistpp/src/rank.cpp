#include <nistpp/tests.h>
#include <nistpp/math_helpers.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include <cmath>

namespace nistpp
{

constexpr double        threshold     = 0.01;
constexpr std::size_t   minimumBits   = 32*32*38;

template<std::size_t M, std::size_t Q>
class BitsMatrix
{
public:
    using matrix_t = Eigen::Matrix<float, M, Q>;
    BitsMatrix(const BitsStorage& data, std::size_t k)
    {
        for(size_t i = 0; i < M; ++i)
        {
            for(size_t j = 0; j < Q; ++j)
            {
                matrix_(i, j) = data[k * (M*Q) + j + i * M];
            }
        }
    }

    std::size_t rank()
    {
        Eigen::FullPivLU<matrix_t> tmp(matrix_);
        return tmp.rank();
    }

private:
    matrix_t matrix_;
};

constexpr double calcProduct(int32_t r)
{
    double res = 1;
    for(int32_t i = 0; i <= r-1; ++i)
        res *= nistpp::const_pow(1.0 - static_cast<double>(nistpp::const_pow(2, i-32)), 2) / (1.0 - nistpp::const_pow(2, i-r));
    return res;
}

constexpr double pNumber(int32_t r)
{
    return nistpp::const_pow(2, r*(32+32-r) - 32*32) * calcProduct(r);
}

return_t RankTest(const BitsStorage &data)
{
    constexpr double p_32 = pNumber(32);
    constexpr double p_31 = pNumber(31);
    constexpr double p_30 = 1 - (p_32 + p_31);

    const auto numberOfBits = data.NumberOfBits();
    if(numberOfBits < minimumBits)
    {
        throw std::logic_error(std::to_string(minimumBits) + " bit required! In storage: " + std::to_string(data.NumberOfBits()));
    }

    std::size_t N = numberOfBits / (32*32);

    std::size_t F_32 = 0;
    std::size_t F_31 = 0;
    for(std::size_t k = 0; k < N; ++k)
    {
        BitsMatrix<32, 32> matrix(data, k);
        auto rank = matrix.rank();

        if(rank == 32)
        {
            ++F_32;
        }

        if(rank == 31)
        {
            ++F_31;
        }
    }

    std::size_t F_30 = N - (F_31 + F_32);

    double chi_squared = std::pow(F_32 - N*p_32, 2.0)/(N*p_32) +
                         std::pow(F_31 - N*p_31, 2.0)/(N*p_31) +
                         std::pow(F_30 - N*p_30, 2.0)/(N*p_30);

    double P = std::exp(-chi_squared/2.0);

    return {P >= threshold, P};
}

} // namespace nistpp
