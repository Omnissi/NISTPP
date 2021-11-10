#ifndef TEMPLATES_H
#define TEMPLATES_H

#include <cstddef>
#include <nistpp/bits_storage.h>
#include <utility>

namespace nistpp
{

using iterator_t = const uint8_t*;

std::size_t GetNumberOfRows(std::size_t m);

std::pair<iterator_t, iterator_t> GetTemplatesSequence(std::size_t i, std::size_t m);

} // namespace nistpp

#endif // TEMPLATES_H
