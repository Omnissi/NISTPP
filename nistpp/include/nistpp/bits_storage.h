#ifndef BITS_STORAGE_H
#define BITS_STORAGE_H

#include <nistpp/types.h>

#include <vector>
#include <bitset>

namespace nistpp
{

/// @brief Class-helper for Nist tests.
class BitsStorage
{
public:
    using word_t = std::bitset<8>;
    using bits_t = std::vector<word_t>;

    BitsStorage() = default;
    ~BitsStorage() = default;

    /// @brief Class constructer.
    /// @param[in] data Data for tests. (1 byte = 1 byte data)
    BitsStorage(const sequence_t& data);

    /// @brief Set data for tests
    /// @param[in] data Data for tests. (1 byte = 1 bit data)
    /// @warning 1 byte of input sequence = 1 bit data
    void SetSequenceByteBit(const sequence_t& data);

    uint8_t operator[](size_t index) const;

    const bits_t& GetBits() const;

    size_t NumberOfBits() const;

private:
    bits_t bits_;
};

} // namespace nistpp

#endif // BITS_STORAGE_H
