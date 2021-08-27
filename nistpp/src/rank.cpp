#include <nistpp/tests.h>
#include <nistpp/math_helpers.h>

#include <cmath>

namespace nistpp
{

constexpr double        threshold     = 0.01;
constexpr std::size_t   minimumBits   = 32*32*38;

template<std::size_t M, std::size_t Q>
class BitsMatrix
{
public:
    using rows_t    = std::bitset<Q>;
    using matrix_t  = std::array<rows_t, M>;
    BitsMatrix(const BitsStorage& data, std::size_t k)
    {
        const auto ind_k = k * (M*Q);
        for(size_t i = 0; i < M; ++i)
        {
            const auto ind_i = i * Q;
            auto& bits = matrix_[i];
            for(size_t j = 0; j < Q; ++j)
            {
                bits[j] = data[ind_k + j + ind_i];
            }
        }
    }

    rows_t& operator[](std::size_t ind)
    {
        return matrix_[ind];
    }

    void swap_rows(std::size_t a, std::size_t b)
    {
        std::swap(matrix_[a], matrix_[b]);
    }

private:
    matrix_t matrix_;
};

template<std::size_t M, std::size_t Q>
void perform_elementary_row_operations_forward(const std::size_t& i, BitsMatrix<M, Q>& matrix)
{
    auto& bits_i = matrix[i];
    for (std::size_t j = i+1; j < M;  ++j)
    {
        if ( matrix[j][i] == 1 )
        {
            auto& bits = matrix[j];
            for (std::size_t k=i; k < Q; ++k)
            {
                bits[k] = (bits[k] + bits_i[k]) % 2;
            }
        }
    }
}

template<std::size_t M, std::size_t Q>
void perform_elementary_row_operations_backward(const std::size_t& i, BitsMatrix<M, Q>& matrix)
{
    auto& bits_i = matrix[i];
    for(int j = i-1; j >= 0; --j)
    {
        if(matrix[j][i] == 1)
        {
            auto& bits = matrix[j];
            for(std::size_t k = 0; k < Q; ++k)
            {
                bits[k] = (bits[k] + bits_i[k]) % 2;
            }
        }
    }
}

template<std::size_t M, std::size_t Q>
std::size_t find_unit_element_and_swap_forward(const std::size_t& i, BitsMatrix<M, Q>& matrix)
{
    std::size_t row_op = 0;

    std::size_t index = i + 1;
    while((index < M) && (matrix[index][i] == 0))
    {
        ++index;
    }
    if(index < M)
    {
        matrix.swap_rows(i, index);
        return 1;
    }

    return 0;
}

template<std::size_t M, std::size_t Q>
std::size_t find_unit_element_and_swap_backward(const std::size_t& i, BitsMatrix<M, Q>& matrix)
{
    int index = i-1;
    while ((index >= 0) && (matrix[index][i] == 0))
    {
        --index;
    }

    if ( index >= 0 )
    {
        matrix.swap_rows(i, index);
        return 1;
    }

    return 0;
}

template<std::size_t M, std::size_t Q>
std::size_t determine_rank(BitsMatrix<M, Q>& matrix)
{
    /* DETERMINE RANK, THAT IS, COUNT THE NUMBER OF NONZERO ROWS */

    std::size_t rank = std::min(M, Q);
    for(std::size_t i = 0; i < M; ++i)
    {
        bool allZeroes = true;
        auto& bits = matrix[i];
        for (std::size_t j = 0; j < Q; ++j)
        {
            if(bits[j] == 1)
            {
                allZeroes = false;
                break;
            }
        }

        if(allZeroes)
        {
            --rank;
        }
    }

    return rank;
}

template<std::size_t M, std::size_t Q>
std::size_t computeRank(BitsMatrix<M, Q>& matrix)
{
    std::size_t rank = 0;
    constexpr auto m  = std::min(M,Q);

    /* FORWARD APPLICATION OF ELEMENTARY ROW OPERATIONS */
    for(std::size_t i = 0; i < m-1; ++i)
    {
        if(matrix[i][i])
        {
            perform_elementary_row_operations_forward(i, matrix);
        }
        else
        {
            if(find_unit_element_and_swap_forward(i, matrix) == 1)
            {
                perform_elementary_row_operations_forward(i, matrix);
            }
        }
    }

    /* BACKWARD APPLICATION OF ELEMENTARY ROW OPERATIONS */
    for (size_t i= m-1; i > 0; --i )
    {
        if (matrix[i][i] == 1)
        {
            perform_elementary_row_operations_backward(i, matrix);
        }
        else
        {
            if (find_unit_element_and_swap_backward(i, matrix) == 1)
            {
                perform_elementary_row_operations_backward(i, matrix);
            }
        }
    }

    rank = determine_rank(matrix);

    return rank;
}

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
        auto rank = computeRank(matrix);

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
