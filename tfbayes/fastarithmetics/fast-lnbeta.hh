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

#ifndef __TFBAYES_FASTARITHMETICS_FAST_LNBETA_HH__
#define __TFBAYES_FASTARITHMETICS_FAST_LNBETA_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/array.hpp>

#include <tfbayes/fastarithmetics/fast-lngamma.hh>

template <class T>
typename T::value_type fast_lnbeta(
        const T& alpha)
{
        typename T::value_type sum1 = 0;
        typename T::value_type sum2 = 0;

        for (size_t i = 0; i < alpha.size(); i++) {
                sum1 += alpha[i];
                sum2 += fast_lngamma(alpha[i]);
        }

        return sum2 - fast_lngamma(sum1);
}

template <class T>
typename T::value_type fast_lnbeta(
        const T& counts,
        const T& alpha)
{
        assert(counts.size() == alpha.size());
        typename T::value_type sum1 = 0;
        typename T::value_type sum2 = 0;

        for (size_t i = 0; i < alpha.size(); i++) {
                sum1 += counts[i] + alpha[i];
                sum2 += fast_lngamma(counts[i] + alpha[i]);
        }

        return sum2 - fast_lngamma(sum1);
}

#endif /* __TFBAYES_FASTARITHMETICS_FAST_LNBETA_HH__ */
