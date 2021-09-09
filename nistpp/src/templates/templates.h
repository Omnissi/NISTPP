#ifndef TEMPLATES_H
#define TEMPLATES_H

#include <nistpp/bits_storage.h>

#include <map>
#include <functional>

namespace nistpp
{

using templ_func_t = std::function<bool(std::size_t, BitsStorage::bits_t::const_iterator, BitsStorage::bits_t::const_iterator)>;

templ_func_t GetTemplatesFunction(std::size_t m, std::size_t &numberOfRows);

} // namespace nistpp

#endif // TEMPLATES_H
