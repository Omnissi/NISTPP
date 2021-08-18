#include <nistpp/bits_storage.h>

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
    for(auto it = tmp.begin(); it < tmp.end(); tmp += 8)
    {
        bits_.emplace_back(std::string(it, it + 8));
    }
}

uint8_t BitsStorage::operator[](size_t index) const
{
    auto div = std::div(index, 8);

    return bits_[div.quot][div.rem];
}

const BitsStorage::bits_t &BitsStorage::GetBits() const
{
    return bits_;
}

size_t BitsStorage::NumberOfBits() const
{
    return bits_.size() * 8;
}

} // namespace nistpp
