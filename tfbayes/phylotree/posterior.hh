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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/array.hpp>

#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/utility/statistics.hh>

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
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                result[i] = 0;
        }
        /* posterior expectation */
        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                double sum = 0;
                double tmp;
                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        sum += it->exponent()[i] + alpha[i];
                }
                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        tmp  = log(it->coefficient());
                        tmp += mbeta_log(it->exponent(), alpha);
                        tmp -= norm;
                        tmp  = exp(tmp)*(it->exponent()[i] + alpha[i])/sum;
                        result[i] += tmp;
                }
        }

        return result;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
boost::array<double, ALPHABET_SIZE> pt_posterior_expectation_prime(
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& polynomial,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha)
{
        double norm = -HUGE_VAL;
        boost::array<double, ALPHABET_SIZE> result;

        /* compute normalization constant */
        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                norm = logadd(norm, log(it->coefficient()) +
                              mbeta_log(it->exponent(), alpha,               0, ALPHABET_SIZE/2) +
                              mbeta_log(it->exponent(), alpha, ALPHABET_SIZE/2, ALPHABET_SIZE  ));
        }
        /* initialize result */
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                result[i] = 0;
        }
        /* posterior expectation */
        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                double sum1 = 0;
                double sum2 = 0;
                double tmp;
                for (size_t i = 0; i < ALPHABET_SIZE/2; i++) {
                        sum1 += it->exponent()[i] + alpha[i];
                }
                for (size_t i = ALPHABET_SIZE/2; i < ALPHABET_SIZE; i++) {
                        sum2 += it->exponent()[i] + alpha[i];
                }
                for (size_t i = 0; i < ALPHABET_SIZE/2; i++) {
                        tmp  = log(it->coefficient());
                        tmp += mbeta_log(it->exponent(), alpha,               0, ALPHABET_SIZE/2);
                        tmp += mbeta_log(it->exponent(), alpha, ALPHABET_SIZE/2, ALPHABET_SIZE);
                        tmp -= norm;
                        tmp  = exp(tmp);
                        tmp *= (it->exponent()[              i] + alpha[              i])/sum1;
                        tmp *= (it->exponent()[ALPHABET_SIZE-1] + alpha[ALPHABET_SIZE-1])/sum2;
                        result[i] += tmp;
                }
        }

        return result;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
double pt_posterior_density(
        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& likelihood,
        const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
        const boost::array<double, ALPHABET_SIZE>& theta)
{
        exponent_t<CODE_TYPE, ALPHABET_SIZE> prior;

        /* for the density we need to substract one from alpha */
        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                prior[i] = alpha[i] - 1;
        }

        double result = likelihood.log_eval(theta) + prior.log_eval(theta) - mbeta_log(alpha);

        return exp(result - pt_marginal_likelihood(likelihood, alpha));
}

#endif /* POSTERIOR_HH */
