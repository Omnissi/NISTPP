#include <nistpp/bits_storage.h>

#include <algorithm>

namespace nistpp
{

BitsStorage::BitsStorage(const sequence_t &data)
{
    bits_.reserve(data.size() * numberOfBitsInWord);

    for(const auto& el : data)
    {
        word_t tmp(el);
        ones_ += tmp.count();
        for(ssize_t i = numberOfBitsInWord - 1; i >= 0; --i)
        {
            bits_.emplace_back(tmp[static_cast<size_t>(i)]);
        }
    }
}

void BitsStorage::SetSequenceByteBit(const sequence_t &data)
{
    bits_.clear();
    bits_.reserve(data.size());

    std::copy(data.begin(), data.end(), bits_.begin());

    ones_ = static_cast<size_t>(std::count(bits_.begin(), bits_.end(), true));
}

bool BitsStorage::operator[](std::size_t index) const
{
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
