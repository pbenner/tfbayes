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

#ifndef MARGINAL_LIKELIHOOD_HH
#define MARGINAL_LIKELIHOOD_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/polynomial.hh>

#include <phylotree.hh>
#include <phylotree-polynomial.hh>

#include <math.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double mbeta_log(
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& exponent,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                sum1 += exponent[i] + alpha[i];
                sum2 += gsl_sf_lngamma(exponent[i] + alpha[i]);
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double pt_marginal_likelihood(
        pt_root_t* node,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha)
{
        double result = 0.0;

        const pt_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(node);

        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                result += log(it->coefficient()) + mbeta_log(it->exponent(), alpha);
        }

        return result;
}

#endif /* MARGINAL_LIKELIHOOD_HH */
