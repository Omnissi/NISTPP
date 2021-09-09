#ifndef TEMPLATE_6_H
#define TEMPLATE_6_H
#include <sprout/valarray.hpp>

namespace nistpp
{
static constexpr sprout::valarray<sprout::valarray<uint8_t, 6>, 20> template6=
{
    {0,0,0,0,0,1},
    {0,0,0,0,1,1},
    {0,0,0,1,0,1},
    {0,0,0,1,1,1},
    {0,0,1,0,1,1},
    {0,0,1,1,0,1},
    {0,0,1,1,1,1},
    {0,1,0,0,1,1},
    {0,1,0,1,1,1},
    {0,1,1,1,1,1},
    {1,0,0,0,0,0},
    {1,0,1,0,0,0},
    {1,0,1,1,0,0},
    {1,1,0,0,0,0},
    {1,1,0,0,1,0},
    {1,1,0,1,0,0},
    {1,1,1,0,0,0},
    {1,1,1,0,1,0},
    {1,1,1,1,0,0},
    {1,1,1,1,1,0},
};
} // namespace nistpp

#endif // TEMPLATE_6_H
