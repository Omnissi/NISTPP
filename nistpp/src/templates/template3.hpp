#ifndef TEMPLATE_3_H
#define TEMPLATE_3_H

#include <array>

namespace nistpp
{
static constexpr std::array<std::array<uint8_t, 3>, 4> template3=
{{
    {0,0,1},
    {0,1,1},
    {1,0,0},
    {1,1,0},
}};
} // namespace nistpp

#endif // TEMPLATE_3_H
