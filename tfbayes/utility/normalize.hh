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

#ifndef __TFBAYES_UTILITY_NORMALIZE_HH__
#define __TFBAYES_UTILITY_NORMALIZE_HH__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <numeric> /* accumulate */

#include <boost/foreach.hpp>

#include <tfbayes/utility/logarithmetic.hh>

template <typename T>
T normalize(const T& container)
{
        T tmp(container);
        typename T::value_type sum = std::accumulate(tmp.begin(), tmp.end(), static_cast<typename T::value_type>(0));

        BOOST_FOREACH(typename T::value_type& n, tmp) {
                n /= sum;
        }
        return tmp;
}

template <typename T>
T log_normalize(const T& container)
{
        T tmp(container);
        typename T::value_type sum = -std::numeric_limits<typename T::value_type>::infinity();

        BOOST_FOREACH(typename T::value_type& n, tmp) {
                sum = logadd(sum, n);
        }
        BOOST_FOREACH(typename T::value_type& n, tmp) {
                n -= sum;
        }
        return tmp;
}

#endif /* __TFBAYES_UTILITY_NORMALIZE_HH__ */
