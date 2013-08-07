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

#ifndef STATISTICS_HH
#define STATISTICS_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cstdlib>
#include <cmath>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include <tfbayes/utility/polynomial.hh>

using namespace std;

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double mbeta_log(
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < ALPHABET_SIZE && alpha[i] > 0; i++) {
                sum1 += alpha[i];
                sum2 += gsl_sf_lngamma(alpha[i]);
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double mbeta_log(
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& exponent,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                if (exponent[i] + alpha[i] > 0) {
                        sum1 += exponent[i] + alpha[i];
                        sum2 += gsl_sf_lngamma(exponent[i] + alpha[i]);
                }
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double mbeta_log(
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& exponent,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
        const size_t from, const size_t to)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = from; i < to && alpha[i] > 0; i++) {
                if (exponent[i] + alpha[i] > 0) {
                        sum1 += exponent[i] + alpha[i];
                        sum2 += gsl_sf_lngamma(exponent[i] + alpha[i]);
                }
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
size_t categorial_sample(const std::vector<double>& p)
{
        double r = (double)rand()/RAND_MAX;
        double c = 0.0;

        for (size_t i = 0; i < ALPHABET_SIZE-1; i++) {
                /* comute cumulative probabilities */
                c += p[i];
                /* return if c is large enough */
                if (r <= c) return i;
        }
        return ALPHABET_SIZE-1;
}

template <size_t ALPHABET_SIZE>
std::vector<double> dirichlet_sample(const std::vector<double>& _alpha, gsl_rng* r)
{
        double alpha[ALPHABET_SIZE];
        double theta[ALPHABET_SIZE];
        std::vector<double> _theta(ALPHABET_SIZE, 0);

        /* copy pseudo counts */
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                alpha[i] = _alpha[i];
        }
        /* generate sample */
        gsl_ran_dirichlet(r, ALPHABET_SIZE, alpha, theta);
        /* copy result */
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                _theta[i] = theta[i];
        }
        return _theta;
}

#endif /* STATISTICS_HH */
