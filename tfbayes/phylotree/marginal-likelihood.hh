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

#ifndef __TFBAYES_PHYLOTREE_MARGINAL_LIKELIHOOD_HH__
#define __TFBAYES_PHYLOTREE_MARGINAL_LIKELIHOOD_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <limits>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-polynomial.hh>
#include <tfbayes/phylotree/phylotree-gradient.hh>
#include <tfbayes/utility/statistics.hh>
#include <tfbayes/utility/logarithmetic.hh>
#include <tfbayes/utility/polynomial.hh>
#include <tfbayes/utility/thread-pool.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

// marginal likelihood
////////////////////////////////////////////////////////////////////////////////

template <size_t AS, typename PC>
double pt_marginal_likelihood(
        const polynomial_t<AS, PC>& polynomial,
        const exponent_t<AS, PC>& alpha)
{
        double result = -std::numeric_limits<double>::infinity();
        double mbeta_alpha = mbeta_log(alpha);

        for (typename polynomial_t<AS, PC>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                result = logadd(result, log(it->coefficient()) +
                                mbeta_log(it->exponent(), alpha) - mbeta_alpha);
        }

        return result;
}

template <size_t AS, typename AC, typename PC>
double pt_marginal_likelihood(
        const pt_root_t& tree,
        const std::vector<AC>& observations,
        const exponent_t<AS, PC>& alpha)
{
        const polynomial_t<AS, PC> polynomial =
                pt_likelihood<AS, AC, PC>(tree, observations);

        return pt_marginal_likelihood<AS, PC>(polynomial, alpha);
}

template <size_t AS, typename AC, typename PC>
double pt_marginal_likelihood(
        // the alignment contains n columns of this form
        double n,
        const pt_root_t& tree,
        const std::vector<AC>& observations,
        const exponent_t<AS, PC>& alpha)
{
        return n*pt_marginal_likelihood(tree, observations, alpha);
}

template <size_t AS, typename AC, typename PC>
double pt_marginal_likelihood(
        const pt_root_t& tree,
        const alignment_map_t<AC>& alignment,
        const exponent_t<AS, PC>& alpha,
        thread_pool_t& thread_pool
        ) {
        // type of the marginal likelihood function that is called by
        // the thread pool
        typedef double (*mlf)(double, const pt_root_t&, const std::vector<AC>&, const exponent_t<AS, PC>&);
        double result = 0;
        // results for each column are stored in a vector of
        // futures (with pre allocated capacity)
        future_vector_t<double> futures(alignment.size());
        // current position in the alignment
        size_t i = 0;
        // launch threads to compute the likelihood
        for (typename alignment_map_t<AC>::const_iterator it = alignment.begin();
             it != alignment.end(); it++) {
                boost::function<double ()> f = boost::bind(
                        static_cast<mlf>(&pt_marginal_likelihood<AS, AC, PC>),
                        static_cast<double>(it->second),
                        boost::cref(tree),
                        boost::cref(it->first),
                        boost::cref(alpha));
                futures[i++] = thread_pool.schedule(f);
        }
        for (size_t i = 0; i < futures.size(); i++) {
                result += futures[i].get();
        }
        return result;
}

// derivative of the marginal likelihood
////////////////////////////////////////////////////////////////////////////////

template <size_t AS, typename PC>
pt_marginal_derivative_t
pt_marginal_derivative(
        const pt_polynomial_derivative_t<AS, PC>& polynomial,
        const exponent_t<AS, PC>& alpha)
{
        pt_marginal_derivative_t result(polynomial.d());
        double mbeta_alpha = mbeta_log(alpha);

        for (typename polynomial_t<AS, PC>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                result += it->coefficient()*std::exp(mbeta_log(it->exponent(), alpha) - mbeta_alpha);
        }
        for (size_t i = 0; i < polynomial.d(); i++) {
                for (typename polynomial_t<AS, PC>::const_iterator it = polynomial.derivative()[i].begin();
                     it != polynomial.derivative()[i].end(); it++) {
                        result.derivative()[i] +=
                                it->coefficient()*std::exp(mbeta_log(it->exponent(), alpha) - mbeta_alpha);
                }
        }
        return result;
}

template <size_t AS, typename AC, typename PC>
pt_marginal_derivative_t
pt_marginal_derivative(
        const pt_root_t& tree,
        const std::vector<AC>& observations,
        const exponent_t<AS, PC>& alpha)
{
        const pt_polynomial_derivative_t<AS, PC> polynomial =
                pt_likelihood_derivative<AS, AC, PC>(tree, observations);

        return pt_marginal_derivative<AS, PC>(polynomial, alpha);
}

template <size_t AS, typename AC, typename PC>
pt_marginal_derivative_t
pt_marginal_derivative(
        // the alignment contains n columns of this form
        double n,
        const pt_root_t& tree,
        const std::vector<AC>& observations,
        const exponent_t<AS, PC>& alpha)
{
        pt_marginal_derivative_t result =
                pt_marginal_derivative(tree, observations, alpha);

        for (size_t i = 0; i < result.d(); i++) {
                result.derivative()[i] = n*result.derivative()[i]/result;
        }
        // from now on use log probabilities
        result = n*std::log(result);

        return result;
}

template <size_t AS, typename AC, typename PC>
pt_marginal_derivative_t
pt_marginal_derivative(
        const pt_root_t& tree,
        const alignment_map_t<AC>& alignment,
        const exponent_t<AS, PC>& alpha,
        thread_pool_t& thread_pool
        ) {
        // type of the marginal likelihood function that is called by
        // the thread pool
        typedef pt_marginal_derivative_t (*mlf)(
                double, const pt_root_t&, const std::vector<AC>&, const exponent_t<AS, PC>&);
        // resulting derivative
        pt_marginal_derivative_t result(tree.n_nodes-1);
        // results for each column are stored in a vector of
        // futures (with pre allocated capacity)
        future_vector_t<pt_marginal_derivative_t> futures(alignment.size());
        // current position in the alignment
        size_t i = 0;
        // launch threads to compute the likelihood
        for (typename alignment_map_t<AC>::const_iterator it = alignment.begin();
             it != alignment.end(); it++) {
                boost::function<pt_marginal_derivative_t ()> f = boost::bind(
                        static_cast<mlf>(&pt_marginal_derivative<AS, AC, PC>),
                        static_cast<double>(it->second),
                        boost::cref(tree),
                        boost::cref(it->first),
                        boost::cref(alpha));
                futures[i++] = thread_pool.schedule(f);
        }
        // marginal likelihood
        for (size_t i = 0; i < futures.size(); i++) {
                result += futures[i].get();
        }
        // derivatives
        for (ssize_t i = 0; i < result.d(); i++) {
                for (size_t j = 0; j < futures.size(); j++) {
                        result.derivative()[i] += futures[j].get().derivative()[i];
                }
        }
        return result;
}

#endif /* __TFBAYES_PHYLOTREE_MARGINAL_LIKELIHOOD_HH__ */
