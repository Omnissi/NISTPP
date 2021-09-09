#include <nistpp/bits_storage.h>

#include <numeric>

namespace nistpp
{

BitsStorage::BitsStorage(const sequence_t &data)
{
    bits_.reserve(data.size() * numberOfBitsInWord);

    for(const auto& el : data)
    {
        word_t tmp(el);
        ones_ += tmp.count();
        for(int32_t i = numberOfBitsInWord - 1; i >= 0; --i)
        {
            bits_.emplace_back(tmp[i]);
        }
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
        word_t tmp(std::string(it, it + numberOfBitsInWord));
        ones_ += tmp.count();
        for(int32_t i = numberOfBitsInWord - 1; i >= 0; --i)
        {
            bits_.emplace_back(tmp[i]);
        }
    }
}

uint8_t BitsStorage::operator[](size_t index) const
{
//    auto div = std::div(index, static_cast<int32_t>(numberOfBitsInWord));

//    return bits_[div.quot][numberOfBitsInWord - 1 - div.rem];
    return bits_[index];
}

const BitsStorage::bits_t &BitsStorage::GetBits() const
{
    return bits_;
}

size_t BitsStorage::NumberOfBits() const
{
    return bits_.size();
}

size_t BitsStorage::NumberOfOnes() const
{
    return ones_;
}

} // namespace nistpp
