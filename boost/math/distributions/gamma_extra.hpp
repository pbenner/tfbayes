/* Copyright (C) 2012 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef __BOOST_MATH_DISTRIBUTIONS_GAMMA_EXTRA_HPP__
#define __BOOST_MATH_DISTRIBUTIONS_GAMMA_EXTRA_HPP__

#include <boost/math/distributions/gamma.hpp>

namespace boost { namespace math {

template <class RealType, class Policy>
inline RealType pdf_derivative(const gamma_distribution<RealType, Policy>& dist, const RealType& x)
{
        RealType a = dist.shape();
        RealType b = dist.scale();

        return ((a-1)/x - 1/b)*pdf(dist, x);
}

template <class RealType, class Policy>
inline RealType log_pdf_derivative(const gamma_distribution<RealType, Policy>& dist, const RealType& x)
{
        RealType a = dist.shape();
        RealType b = dist.scale();

        return ((a-1)/x - 1/b);
}


} // namespace math
} // namespace boost

#endif /* __BOOST_MATH_DISTRIBUTIONS_GAMMA_EXTRA_HPP__ */
