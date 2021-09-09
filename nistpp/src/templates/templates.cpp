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

namespace nistpp
{

template<class templ_t, class It>
bool templateEqual(std::size_t i, const templ_t templ, It begin, It end)
{
    if(static_cast<size_t>(std::distance(begin, end)) != templ[0].size())
    {
        throw std::logic_error("Input iterators incorrect!");
    }

    return std::equal(templ[i].begin(), templ[i].end(), begin, end);
}

bool templateEqual2(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template2, begin, end);
}

bool templateEqual3(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template3, begin, end);
}

bool templateEqual4(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template4, begin, end);
}

bool templateEqual5(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template5, begin, end);
}

bool templateEqual6(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template6, begin, end);
}

bool templateEqual7(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template7, begin, end);
}

bool templateEqual8(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template8, begin, end);
}

bool templateEqual9(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template9, begin, end);
}

bool templateEqual10(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template10, begin, end);
}

bool templateEqual11(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template11, begin, end);
}

bool templateEqual12(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template12, begin, end);
}

bool templateEqual13(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template13, begin, end);
}

bool templateEqual14(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template14, begin, end);
}

bool templateEqual15(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template15, begin, end);
}

bool templateEqual16(std::size_t i, BitsStorage::bits_t::const_iterator begin, BitsStorage::bits_t::const_iterator end)
{
    return templateEqual(i, template16, begin, end);
}

templ_func_t GetTemplatesFunction(std::size_t m, std::size_t& numberOfRows)
{
    switch (m)
    {
    case 2: numberOfRows = template2.size(); return std::bind(templateEqual2, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 3: numberOfRows = template3.size(); return std::bind(templateEqual3, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 4: numberOfRows = template4.size(); return std::bind(templateEqual4, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 5: numberOfRows = template5.size(); return std::bind(templateEqual5, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 6: numberOfRows = template6.size(); return std::bind(templateEqual6, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 7: numberOfRows = template7.size(); return std::bind(templateEqual7, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 8: numberOfRows = template8.size(); return std::bind(templateEqual8, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 9: numberOfRows = template9.size(); return std::bind(templateEqual9, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 10: numberOfRows = template10.size(); return std::bind(templateEqual10, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 11: numberOfRows = template11.size(); return std::bind(templateEqual11, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 12: numberOfRows = template12.size(); return std::bind(templateEqual12, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 13: numberOfRows = template13.size(); return std::bind(templateEqual13, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 14: numberOfRows = template14.size(); return std::bind(templateEqual14, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 15: numberOfRows = template15.size(); return std::bind(templateEqual15, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    case 16: numberOfRows = template16.size(); return std::bind(templateEqual16, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    }

    return nullptr;
}

} // namespace nistpp
