#ifndef TEMPLATE_4_H
#define TEMPLATE_4_H
#include <sprout/valarray.hpp>

namespace nistpp
{
static constexpr sprout::valarray<sprout::valarray<bool, 4>, 6> template4=
{
    {0,0,0,1},
    {0,0,1,1},
    {0,1,1,1},
    {1,0,0,0},
    {1,1,0,0},
    {1,1,1,0},
};
} // namespace nistpp

#endif // TEMPLATE_4_H
