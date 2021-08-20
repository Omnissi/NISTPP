#include <nistpp/bits_storage.h>

#include <numeric>

namespace nistpp
{

BitsStorage::BitsStorage(const sequence_t &data)
{
    bits_.reserve(data.size());

    for(const auto& el : data)
    {
        bits_.emplace_back(el);
    }
}

void BitsStorage::SetSequenceByteBit(const sequence_t &data)
{
    //TODO: refactor

    std::string tmp(data.begin(), data.end());

    bits_.clear();
    bits_.reserve(tmp.size());
    for(auto it = tmp.begin(); it < tmp.end(); tmp += numberOfBitsInWord)
    {
        bits_.emplace_back(std::string(it, it + numberOfBitsInWord));
    }
}

uint8_t BitsStorage::operator[](size_t index) const
{
    auto div = std::div(index, static_cast<int32_t>(numberOfBitsInWord));

    return bits_[div.quot][numberOfBitsInWord - 1 - div.rem];
}

const BitsStorage::bits_t &BitsStorage::GetBits() const
{
    return bits_;
}

size_t BitsStorage::NumberOfBits() const
{
    return bits_.size() * 8;
}

size_t BitsStorage::NumberOfOnes() const
{
    static auto ones = std::accumulate(bits_.begin(), bits_.end(), 0,
                                [](const auto& a, const auto& b)
                                {
                                    return a + b.count();
                                });

    return static_cast<size_t>(ones);
}

} // namespace nistpp
