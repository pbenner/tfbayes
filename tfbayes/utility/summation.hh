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

#ifndef __TFBAYES_UTILITY_SUMMATION_HH__
#define __TFBAYES_UTILITY_SUMMATION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <boost/foreach.hpp>

template <class RealType>
RealType GCC_ATTRIBUTE_NOAMATH
ksum(const std::vector<RealType>& input)
{
        // Standard Kahan summation
        RealType sum = 0.0;
        RealType c   = 0.0;
        for (size_t i = 0; i < input.size(); i++) {
                RealType y = input[i] - c;
                RealType t = sum + y;
                c = (t - sum) - y;
                sum = t;
        }
        return sum;
}

template <class RealType>
RealType GCC_ATTRIBUTE_NOAMATH
msum(const std::vector<RealType>& input) {
        // Code from: http://code.activestate.com/recipes/393090/
        // Full precision summation using multiple floats for intermediate values
        // Rounded x+y stored in hi with the round-off stored in lo.  Together
        // hi+lo are exactly equal to x+y.  The inner loop applies hi/lo summation
        // to each partial so that the list of partial sums remains exact.
        // Depends on IEEE-754 arithmetic guarantees.  See proof of correctness at:
        // www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
        RealType hi;
        RealType lo;
        std::vector<RealType> partials;
        std::vector<RealType> partials_tmp;
        partials.reserve(input.size());
        partials_tmp.reserve(input.size());

        BOOST_FOREACH(RealType x, input) {
                size_t i = 0;
                BOOST_FOREACH(const RealType& y, partials) {
                        if (std::abs(x) < std::abs(y)) {
                                hi = x + y;
                                lo = x - (hi - y);
                        }
                        else {
                                hi = x + y;
                                lo = y - (hi - x);
                        }
                        if (lo) {
                                if (i < partials_tmp.size()) {
                                        partials_tmp[i] = lo;
                                }
                                else {
                                        partials_tmp.push_back(lo);
                                }
                                i++;
                        }
                        x = hi;
                }
                partials_tmp.resize(i);
                partials_tmp.push_back(x);
                partials = partials_tmp;
        }
        return std::accumulate(partials.begin(), partials.end(), RealType(0));
}

#endif /* __TFBAYES_UTILITY_SUMMATION_HH__ */
