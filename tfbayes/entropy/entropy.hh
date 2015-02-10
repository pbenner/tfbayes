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

#ifndef __TFBAYES_ENTROPY_ENTROPY_HH__
#define __TFBAYES_ENTROPY_ENTROPY_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cmath>

#include <boost/foreach.hpp>

#include <tfbayes/utility/summation.hh>

template <class RealType>
RealType entropy(const std::vector<RealType>& v) {
        RealType result = 0.0;
        BOOST_FOREACH(const RealType& x, v) {
                result += -x*static_cast<RealType>(std::log(x));
        }
        return result;
}

template <class RealType>
RealType entropy(const std::vector<RealType>& v, bool extended_precision) {
        if (!extended_precision) {
                return entropy<RealType>(v);
        }
        std::vector<RealType> result(v.size(), 0.0);
        for (size_t i = 0; i < v.size(); i++) {
                result[i] = -v[i]*static_cast<RealType>(std::log(v[i]));
        }
        return msum(result);
}

#endif /* __TFBAYES_ENTROPY_ENTROPY_HH__ */
