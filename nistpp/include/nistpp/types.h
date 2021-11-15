#ifndef NISTPP_TYPES_H
#define NISTPP_TYPES_H

#pragma once

#include <cstdint>
#include <vector>
#include <utility>

namespace nistpp
{

/// @brief Return type all tests.
/// @details First:  true - test passed, else - false
///          Second: P-value
using return_t = std::pair<bool, double>;

/// @brief Type of sequence.
using sequence_t = std::vector<uint8_t>;

/// @brief If P-value < threshold then sequence is not random. Equal ALPHA.
constexpr double threshold = 0.01;

} // namespace nistpp

#endif // NISTPP_TYPES_H
