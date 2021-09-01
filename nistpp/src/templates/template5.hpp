#ifndef TEMPLATE_5_H
#define TEMPLATE_5_H
#include <sprout/valarray.hpp>

namespace nistpp
{
static constexpr sprout::valarray<sprout::valarray<bool, 5>, 12> template5=
{
    {0,0,0,0,1},
    {0,0,0,1,1},
    {0,0,1,0,1},
    {0,0,1,1,1},
    {0,1,0,1,1},
    {0,1,1,1,1},
    {1,0,0,0,0},
    {1,0,1,0,0},
    {1,1,0,0,0},
    {1,1,0,1,0},
    {1,1,1,0,0},
    {1,1,1,1,0},
};
} // namespace nistpp

#endif // TEMPLATE_5_H
