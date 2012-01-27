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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <likelihood.hh>

#include <iostream>
#include <gsl/gsl_sf_gamma.h>

using namespace std;

double mbeta_log(
        const polynomial_term_t<code_t, alphabet_size>& term,
        const polynomial_term_t<code_t, alphabet_size>& alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < alphabet_size; i++) {
                sum1 += term[i] + alpha[i];
                sum2 += gsl_sf_lngamma(term[i] + alpha[i]);
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

double
pt_marginal_likelihood(
        pt_root_t<code_t, alphabet_size>* node,
        const polynomial_term_t<code_t, alphabet_size> alpha)
{
        double result = 0.0;

        pt_likelihood(node);
        const polynomial_t<code_t, alphabet_size>& polynomial = node->poly_sum;

        for (polynomial_t<code_t, alphabet_size>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                result += log(it->coefficient()) + mbeta_log(*it, alpha);
        }

        return result;
}
