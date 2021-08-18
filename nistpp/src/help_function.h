#ifndef HELP_FUNCTION_H
#define HELP_FUNCTION_H

#include <cstddef>
#include <utility>

namespace nistpp
{
/// @brief Calculating 1 in word.
/// @param[in] x Word.
/// @return Number of 1 in word.
/// @tparam T Type of word.
template<class T>
inline std::size_t CalcOnes(T x)
{
    static_assert ( std::is_integral<T>::value, "");

    constexpr std::size_t   numOfBytes    = sizeof (T);
    constexpr std::size_t   degree        = (numOfBytes == 1) ? 0    : (numOfBytes == 2 ? 1      : (numOfBytes == 4 ? 2 : 3));
    constexpr T             firstMask     = (numOfBytes == 1) ? 0x55 : (numOfBytes == 2 ? 0x5555 : (numOfBytes == 4 ? 0x55555555 : 0x5555555555555555));
    constexpr T             secondMask    = (numOfBytes == 1) ? 0x33 : (numOfBytes == 2 ? 0x3333 : (numOfBytes == 4 ? 0x33333333 : 0x3333333333333333));
    constexpr T             thirdMask     = (numOfBytes == 1) ? 0x0f : (numOfBytes == 2 ? 0x0f0f : (numOfBytes == 4 ? 0x0f0f0f0f : 0x0f0f0f0f0f0f0f0f));

    x = x - ((x >> 1) & firstMask);
    x = (x & secondMask) + ((x >> 2) & secondMask);
    x = (x + (x >> 4)) & thirdMask;

    std::size_t bitsMove = 8;
    for(std::size_t i = 0; i < degree; ++i)
    {
        x = x + ( x >> bitsMove );
        bitsMove <<= 1;
    }

    return x & 0x3f;
}

} // namespace nistpp

#endif // HELP_FUNCTION_H
