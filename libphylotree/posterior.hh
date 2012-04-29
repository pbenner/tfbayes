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

#ifndef POSTERIOR_HH
#define POSTERIOR_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/array.hpp>

#include <marginal-likelihood.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
boost::array<double, ALPHABET_SIZE> pt_posterior_expectation(
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& polynomial,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha)
{
        double norm = -HUGE_VAL;
        boost::array<double, ALPHABET_SIZE> result;

        /* compute normalization constant */
        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                norm = logadd(norm, log(it->coefficient()) + mbeta_log(it->exponent(), alpha));
        }
        /* initialize result */
        for (CODE_TYPE i = 0; i < ALPHABET_SIZE; i++) {
                result[i] = 0;
        }
        /* posterior expectation */
        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                double sum = 0;
                double tmp;
                for (CODE_TYPE i = 0; i < ALPHABET_SIZE; i++) {
                        sum += it->exponent()[i] + alpha[i];
                }
                for (CODE_TYPE i = 0; i < ALPHABET_SIZE; i++) {
                        tmp  = log(it->coefficient());
                        tmp += mbeta_log(it->exponent(), alpha);
                        tmp -= norm;
                        tmp  = exp(tmp)*(it->exponent()[i] + alpha[i])/sum;
                        result[i] += tmp;
                }
        }

        return result;
}


#endif /* POSTERIOR_HH */