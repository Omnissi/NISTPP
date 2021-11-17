#include <nistpp/math_helpers.h>

#include <boost/math/special_functions/gamma.hpp>

namespace nistpp
{

double igamc(double a, double z)
{
    return  boost::math::gamma_q(a, z);
}

} // namespace nistpp