#ifndef TEMPLATE_2_H
#define TEMPLATE_2_H
#include <sprout/valarray.hpp>

namespace nistpp
{
static constexpr sprout::valarray<sprout::valarray<uint8_t, 2>, 2> template2=
{
    {0,1},
    {1,0},
};
} // namespace nistpp

#endif // TEMPLATE_2_H
