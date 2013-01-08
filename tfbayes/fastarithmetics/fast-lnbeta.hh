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

#ifndef FAST_LNBETA_HH
#define FAST_LNBETA_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/array.hpp>

#include <tfbayes/fastarithmetics/fast-lngamma.hh>

template <size_t SIZE>
double fast_lnbeta(
        const boost::array<double, SIZE> alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < SIZE; i++) {
                sum1 += alpha[i];
                sum2 += fast_lngamma(alpha[i]);
        }

        return sum2 - fast_lngamma(sum1);
}

template <size_t SIZE>
double fast_lnbeta(
        const boost::array<double, SIZE> extra,
        const boost::array<double, SIZE> alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < SIZE; i++) {
                sum1 += extra[i] + alpha[i];
                sum2 += fast_lngamma(extra[i] + alpha[i]);
        }

        return sum2 - fast_lngamma(sum1);
}

template <size_t SIZE>
double fast_lnbeta(
        const double* alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < SIZE; i++) {
                sum1 += alpha[i];
                sum2 += fast_lngamma(alpha[i]);
        }

        return sum2 - fast_lngamma(sum1);
}

template <size_t SIZE>
double fast_lnbeta(
        const double* extra,
        const double* alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < SIZE; i++) {
                sum1 += extra[i] + alpha[i];
                sum2 += fast_lngamma(extra[i] + alpha[i]);
        }

        return sum2 - fast_lngamma(sum1);
}

template <size_t SIZE>
double fast_lnbeta(
        const boost::array<double, SIZE> extra,
        const double* alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < SIZE; i++) {
                sum1 += extra[i] + alpha[i];
                sum2 += fast_lngamma(extra[i] + alpha[i]);
        }

        return sum2 - fast_lngamma(sum1);
}

#endif /* FAST_LNBETA_HH */
