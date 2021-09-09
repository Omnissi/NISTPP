#ifndef TEMPLATE_3_H
#define TEMPLATE_3_H
#include <sprout/valarray.hpp>

namespace nistpp
{
static constexpr sprout::valarray<sprout::valarray<uint8_t, 3>, 4> template3=
{
    {0,0,1},
    {0,1,1},
    {1,0,0},
    {1,1,0},
};
} // namespace nistpp

#endif // TEMPLATE_3_H
