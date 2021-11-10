#ifndef TEMPLATE_5_H
#define TEMPLATE_5_H

#include <array>

namespace nistpp
{
static constexpr std::array<std::array<uint8_t, 5>, 12> template5=
{{
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
}};
} // namespace nistpp

#endif // TEMPLATE_5_H
