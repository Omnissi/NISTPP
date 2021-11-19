#ifndef NISTPPP_BITS_STORAGE_H
#define NISTPPP_BITS_STORAGE_H

#pragma once

#include <nistpp/types.h>

#include <vector>
#include <bitset>
#include <memory>

namespace nistpp
{

/// @brief Class-helper for Nist tests.
class BitsStorage
{
public:
    static constexpr std::size_t numberOfBitsInWord = 8;
    using word_t = std::bitset<numberOfBitsInWord>;
    using bits_t = std::vector<bool>;

    BitsStorage() = default;
    ~BitsStorage() = default;

    BitsStorage(const BitsStorage&) = delete;
    BitsStorage(BitsStorage&&) noexcept = delete;

    /// @brief Class constructer.
    /// @param[in] data Data for tests. (1 byte = 1 byte data)
    explicit BitsStorage(const sequence_t& data);

    /// @brief Set data for tests
    /// @param[in] data Data for tests. (1 byte = 1 bit data)
    /// @warning 1 byte of input sequence = 1 bit data
    void SetSequenceByteBit(const sequence_t& data);

    /// @brief Get bit from index.
    /// @param[in] index Index of bit.
    /// @return Bit.
    bool operator[](std::size_t index) const;

    /// @brief Get raw sequence, where 1 byte = 1 bit.
    /// @return Raw sequence.
    const bits_t& GetBits() const;

    /// @brief Get number of bit in sequence.
    /// @return Number of bit in sequence.
    std::size_t NumberOfBits() const;

    /// @brief Get number of 1-s in sequence.
    /// @return Number of 1-s in sequence.
    std::size_t NumberOfOnes() const;

private:
    bits_t bits_{};
    std::size_t ones_{};
};

} // namespace nistpp

#endif // NISTPPP_BITS_STORAGE_H
