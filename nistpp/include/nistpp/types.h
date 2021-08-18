#ifndef TYPES_H
#define TYPES_H

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

} // namespace nistpp

#endif // TYPES_H
