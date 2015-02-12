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

#ifndef __TFBAYES_UTILITY_MULTINOMIAL_BETA_HH__
#define __TFBAYES_UTILITY_MULTINOMIAL_BETA_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cmath>

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

template <class T>
typename T::value_type mbeta_log(
        const T& alpha)
{
        typename T::value_type sum1 = 0;
        typename T::value_type sum2 = 0;

        for (size_t i = 0; i < alpha.size(); i++) {
                sum1 += alpha[i];
                sum2 += boost::math::lgamma<typename T::value_type>(alpha[i]);
        }

        return sum2 - boost::math::lgamma<typename T::value_type>(sum1);
}

template <class T>
typename T::value_type mbeta_log(
        const T& counts,
        const T& alpha)
{
        assert(counts.size() == alpha.size());
        typename T::value_type sum1 = 0;
        typename T::value_type sum2 = 0;

        for (size_t i = 0; i < counts.size(); i++) {
                sum1 += counts[i] + alpha[i];
                sum2 += boost::math::lgamma<typename T::value_type>(counts[i] + alpha[i]);
        }

        return sum2 - boost::math::lgamma<typename T::value_type>(sum1);
}

template <class T>
typename T::value_type mbeta_log(
        const T& counts,
        const T& alpha,
        const size_t from, const size_t to)
{
        assert(counts.size() == alpha.size());
        typename T::value_type sum1 = 0;
        typename T::value_type sum2 = 0;

        for (size_t i = from; i < to && alpha[i] > 0; i++) {
                if (counts[i] + alpha[i] > 0) {
                        sum1 += counts[i] + alpha[i];
                        sum2 += boost::math::lgamma<typename T::value_type>(counts[i] + alpha[i]);
                }
        }

        return sum2 - boost::math::lgamma<typename T::value_type>(sum1);
}

#endif /* __TFBAYES_UTILITY_MULTINOMIAL_BETA_HH__ */
