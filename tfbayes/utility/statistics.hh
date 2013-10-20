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

template <size_t AS, typename PC>
double mbeta_log(
        const exponent_t<AS, PC>& alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < AS && alpha[i] > 0; i++) {
                sum1 += alpha[i];
                sum2 += gsl_sf_lngamma(alpha[i]);
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

template <size_t AS, typename PC>
double mbeta_log(
        const exponent_t<AS, PC>& exponent,
        const exponent_t<AS, PC>& alpha)
{
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i = 0; i < AS; i++) {
                if (exponent[i] + alpha[i] > 0) {
                        sum1 += exponent[i] + alpha[i];
                        sum2 += gsl_sf_lngamma(exponent[i] + alpha[i]);
                }
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

template <size_t AS, typename PC>
double mbeta_log(
        const exponent_t<AS, PC>& exponent,
        const exponent_t<AS, PC>& alpha,
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

template <size_t AS, typename PC>
size_t categorial_sample(const std::vector<double>& p)
{
        double r = (double)rand()/RAND_MAX;
        double c = 0.0;

        for (size_t i = 0; i < AS-1; i++) {
                /* comute cumulative probabilities */
                c += p[i];
                /* return if c is large enough */
                if (r <= c) return i;
        }
        return AS-1;
}

template <size_t AS>
std::vector<double> dirichlet_sample(const std::vector<double>& _alpha, gsl_rng* r)
{
        double alpha[AS];
        double theta[AS];
        std::vector<double> _theta(AS, 0);

        /* copy pseudo counts */
        for (size_t i = 0; i < AS; i++) {
                alpha[i] = _alpha[i];
        }
        /* generate sample */
        gsl_ran_dirichlet(r, AS, alpha, theta);
        /* copy result */
        for (size_t i = 0; i < AS; i++) {
                _theta[i] = theta[i];
        }
        return _theta;
}

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

static inline
size_t select_component(size_t k, double log_weights[], boost::random::mt19937& gen)
{
        /* log_weights are cumulative */
        boost::random::uniform_01<> dist;
        
        const double r = dist(gen);

        if (r == 0.0) {
                return 0;
        }

        const double log_r = log(r) + log_weights[k-1];
        for (size_t i = 0; i < k-1; i++) {
                if (log_r <= log_weights[i]) {
                        return i;
                }
        }
        return k-1;
}

static inline
size_t select_max_component(size_t k, double log_weights[])
{
        /* log_weights are cumulative, so we need to normalize by log_weights[k-1] */

        /* resulting cluster number */
        size_t result = 0;
        /* difference between two successive log weights */
        double diff   = 0.0;
        /* we need to store the last value because zero is not
         * represented in the log weights array */
        double last   = 0.0;

        for (size_t i = 0; i < k; i++) {
                /* if the increase for this cluster is higher than any
                 * before */
                if (diff < exp(log_weights[i] - log_weights[k-1]) - last) {
                        /* store the new increase in probability */
                        diff   = exp(log_weights[i] - log_weights[k-1]) - last;
                        /* and the cluster number */
                        result = i;
                }
                /* also save the value of the last weight */
                last = exp(log_weights[i] - log_weights[k-1]);
        }
        return result;
}

#endif /* STATISTICS_HH */