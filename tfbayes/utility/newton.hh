/* Copyright (C) 2015 Philipp Benner
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

#ifndef __TFBAYES_UTILITY_NEWTON_HH__
#define __TFBAYES_UTILITY_NEWTON_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cmath>

#include <boost/function.hpp>

/* Given a function f, compute the inverse x for a given
 * value y.
 */

template <class RealType>
RealType newton(boost::function<RealType (RealType)> f,
                boost::function<RealType (RealType)> df,
                RealType x, const RealType y)
{
        /* x: current position
         * y: target value
         */
        while (std::abs(f(x) - y) > 1e-8) {
                x += (y - f(x))/df(x);
        }
        return x;
}

#endif /* __TFBAYES_UTILITY_NEWTON_HH__ */
