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

#ifndef __TFBAYES_UTILITY_STATISTICS_HH__
#define __TFBAYES_UTILITY_STATISTICS_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cstdlib>
#include <cmath>

#include <gsl/gsl_randist.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <tfbayes/utility/polynomial.hh>

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

template <size_t AS, typename PC>
size_t categorical_sample(const std::vector<double>& p, boost::random::mt19937& gen)
{
        boost::random::uniform_01<> dist;

        double r = dist(gen);
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

static inline
size_t select_component(size_t k, const double log_weights[], boost::random::mt19937& gen)
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
size_t select_component(const std::vector<double>& log_weights, boost::random::mt19937& gen)
{
        /* log_weights are cumulative */
        boost::random::uniform_01<> dist;
        
        const double r = dist(gen);
        const size_t k = log_weights.size();

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
std::pair<size_t,size_t> select_component2(
        size_t k,
        const double log_weights1[],
        const double log_weights2[],
        boost::random::mt19937& gen)
{
        /* log_weights are cumulative */
        boost::random::uniform_01<> dist;
        
        const double r = dist(gen);

        if (r == 0.0) {
                return std::make_pair(1,0);
        }

        const double log_r = log(r) + log_weights2[k-1];
        for (size_t i = 0; i < k; i++) {
                if (log_r <= log_weights1[i]) {
                        return std::make_pair(1, i);
                }
        }
        for (size_t i = 0; i < k-1; i++) {
                if (log_r <= log_weights2[i]) {
                        return std::make_pair(2, i);
                }
        }
        return std::make_pair(2, k-1);
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

static inline
std::pair<size_t, size_t> select_max_component2(size_t k, const double log_weights1[], const double log_weights2[])
{
        /* log_weights are cumulative, so we need to normalize by log_weights[k-1] */

        /* resulting cluster number */
        std::pair<size_t, size_t> result = std::make_pair(1, 0);
        /* difference between two successive log weights */
        double diff   = 0.0;
        /* we need to store the last value because zero is not
         * represented in the log weights array */
        double last   = 0.0;

        for (size_t i = 0; i < k; i++) {
                /* if the increase for this cluster is higher than any
                 * before */
                if (diff < exp(log_weights1[i] - log_weights2[k-1]) - last) {
                        /* store the new increase in probability */
                        diff   = exp(log_weights1[i] - log_weights2[k-1]) - last;
                        /* and the cluster number */
                        result.first  = 1;
                        result.second = i;
                }
                /* also save the value of the last weight */
                last = exp(log_weights1[i] - log_weights2[k-1]);
        }
        for (size_t i = 0; i < k; i++) {
                /* if the increase for this cluster is higher than any
                 * before */
                if (diff < exp(log_weights2[i] - log_weights2[k-1]) - last) {
                        /* store the new increase in probability */
                        diff   = exp(log_weights2[i] - log_weights2[k-1]) - last;
                        /* and the cluster number */
                        result.first  = 2;
                        result.second = i;
                }
                /* also save the value of the last weight */
                last = exp(log_weights2[i] - log_weights2[k-1]);
        }
        return result;
}

#endif /* __TFBAYES_UTILITY_STATISTICS_HH__ */
