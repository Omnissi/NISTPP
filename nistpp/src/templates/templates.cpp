#include "templates.h"

#include "template2.hpp"
#include "template3.hpp"
#include "template4.hpp"
#include "template5.hpp"
#include "template6.hpp"
#include "template7.hpp"
#include "template8.hpp"
#include "template9.hpp"
#include "template10.hpp"
#include "template11.hpp"
#include "template12.hpp"
#include "template13.hpp"
#include "template14.hpp"
#include "template15.hpp"
#include "template16.hpp"
#include <cstddef>

namespace nistpp
{

std::pair<iterator_t, iterator_t> GetTemplatesSequence(std::size_t i, std::size_t m)
{
#define ret_it(templ) { templ[i].begin(), templ[i].end() }
    switch (m)
    {
    case 2: return ret_it(template2);
    case 3: return ret_it(template3);
    case 4: return ret_it(template4);
    case 5: return ret_it(template5);
    case 6: return ret_it(template6);
    case 7: return ret_it(template7);
    case 8: return ret_it(template8);
    case 9: return ret_it(template9);
    case 10: return ret_it(template10);
    case 11: return ret_it(template11);
    case 12: return ret_it(template12);
    case 14: return ret_it(template14);
    case 15: return ret_it(template15);
    case 16: return ret_it(template16);
    }
#undef ret_it

    return {nullptr, nullptr};
}

std::size_t GetNumberOfRows(std::size_t m)
{
    switch (m)
    {
    case 2: return template2.size();
    case 3: return template3.size();
    case 4: return template4.size();
    case 5: return template5.size();
    case 6: return template6.size();
    case 7: return template7.size();
    case 8: return template8.size();
    case 9: return template9.size();
    case 10: return template10.size();
    case 11: return template11.size();
    case 12: return template12.size();
    case 13: return template13.size();
    case 14: return template14.size();
    case 15: return template15.size();
    case 16: return template16.size();
    }

    return 0;
}

} // namespace nistpp
