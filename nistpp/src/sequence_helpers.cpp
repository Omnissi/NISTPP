#include <nistpp/sequence_helpers.h>

#include <stdexcept>
#include <bitset>
namespace nistpp
{

constexpr std::size_t bitsInByte = 8;

sequence_t ConvertStringToSequence(const std::string &data)
{
    if(data.size() % bitsInByte)
    {
        throw std::invalid_argument("Number of bits must % 8");
    }

    sequence_t res;
    res.reserve(data.size() / bitsInByte);
    auto it = data.begin();
    while(it != data.end())
    {
        std::bitset<bitsInByte> bits(std::string(it, it+bitsInByte));
        res.emplace_back(bits.to_ulong());
        it += bitsInByte;
    }

    return res;
}

} // namespace nistpp
