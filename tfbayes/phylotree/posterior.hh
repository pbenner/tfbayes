/* Copyright (C) 2012-2013 Philipp Benner
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

#ifndef __TFBAYES_PHYLOTREE_POSTERIOR_HH__
#define __TFBAYES_PHYLOTREE_POSTERIOR_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <limits>

#include <boost/array.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/marginal-likelihood.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/distribution.hh>
#include <tfbayes/utility/statistics.hh>
#include <tfbayes/utility/thread-pool.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

// posterior distribution of branch lengths given a full alignment
////////////////////////////////////////////////////////////////////////////////

template <size_t AS, typename AC, typename PC>
double pt_posterior(
        const pt_root_t& tree,
        const alignment_map_t<AC>& alignment,
        const exponent_t<AS, PC>& alpha,
        const boost::math::gamma_distribution<>& gamma_prior,
        thread_pool_t& thread_pool
        ) {
        double result = pt_marginal_likelihood(tree, alignment, alpha, thread_pool);
        // prior on branch lengths
        for (pt_node_t::nodes_t::const_iterator it = tree.begin_nodes();
             it != tree.end_nodes(); it++) {
                // skip the root
                if ((*it)->root()) {
                        continue;
                }
                if ((*it)->d > 0) {
                        result += std::log(boost::math::pdf(gamma_prior, (*it)->d));
                }
        }
        return result;
}

template <size_t AS, typename AC, typename PC>
pt_marginal_derivative_t
pt_posterior_derivative(
        const pt_root_t& tree,
        const alignment_map_t<AC>& alignment,
        const exponent_t<AS, PC>& alpha,
        const boost::math::gamma_distribution<>& gamma_prior,
        thread_pool_t& thread_pool
        ) {
        pt_marginal_derivative_t result = pt_marginal_derivative(tree, alignment, alpha, thread_pool);
        // prior on branch lengths
        for (pt_node_t::nodes_t::const_iterator it = tree.begin_nodes();
             it != tree.end_nodes(); it++) {
                // skip the root
                if ((*it)->root()) {
                        continue;
                }
                // posterior value
                if ((*it)->d > 0) {
                        result += std::log(boost::math::pdf(gamma_prior, (*it)->d));
                }
                // posterior derivative
                result.derivative()[(*it)->id]
                        += boost::math::log_pdf_derivative(gamma_prior, (*it)->d);
        }
        return result;
}

// posterior distribution of one column
////////////////////////////////////////////////////////////////////////////////

template <size_t AS, typename PC>
boost::array<double, AS> pt_posterior_expectation(
        const polynomial_t<AS, PC>& polynomial,
        const exponent_t<AS, PC>& alpha)
{
        double norm = -std::numeric_limits<double>::infinity();
        boost::array<double, AS> result;

        /* compute normalization constant */
        for (typename polynomial_t<AS, PC>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                norm = logadd(norm, log(it->coefficient()) + mbeta_log(it->exponent(), alpha));
        }
        /* initialize result */
        for (size_t i = 0; i < AS; i++) {
                result[i] = 0;
        }
        /* posterior expectation */
        for (typename polynomial_t<AS, PC>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                double sum = 0;
                double tmp;
                for (size_t i = 0; i < AS; i++) {
                        sum += it->exponent()[i] + alpha[i];
                }
                for (size_t i = 0; i < AS; i++) {
                        tmp  = log(it->coefficient());
                        tmp += mbeta_log(it->exponent(), alpha);
                        tmp -= norm;
                        tmp  = exp(tmp)*(it->exponent()[i] + alpha[i])/sum;
                        result[i] += tmp;
                }
        }

        return result;
}

template <size_t AS, typename AC, typename PC>
boost::array<double, AS> pt_posterior_expectation(
        const pt_root_t& node,
        const std::vector<AC>& observations,
        const exponent_t<AS, AC>& alpha)
{
        polynomial_t<AS, PC> polynomial = pt_likelihood<AS, AC, PC>(node, observations);

        return pt_posterior_expectation<AS, PC>(polynomial, alpha);
}

template <size_t AS, typename PC>
double pt_posterior_density(
        const polynomial_t<AS, PC>& likelihood,
        const exponent_t<AS, PC>& alpha,
        const boost::array<double, AS>& theta)
{
        exponent_t<AS, PC> prior;

        /* for the density we need to substract one from alpha */
        for (size_t i = 0; i < AS; i++) {
                prior[i] = alpha[i] - 1;
        }

        double result = likelihood.log_eval(theta) + prior.log_eval(theta) - mbeta_log(alpha);

        return exp(result - pt_marginal_likelihood(likelihood, alpha));
}

#endif /* __TFBAYES_PHYLOTREE_POSTERIOR_HH__ */
