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

#ifndef PHYLOTREE_GRADIENT_HH
#define PHYLOTREE_GRADIENT_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <sys/time.h>

#include <set>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <algorithm> /* std::min */

#include <tfbayes/exception.h>
#include <tfbayes/polynomial.hh>

#include <alignment.hh>
#include <phylotree.hh>
#include <phylotree-gradient-coefficient.hh>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>

class gamma_distribution_t {
public:
        gamma_distribution_t(double r, double lambda)
                : r(r), r_gamma(gsl_sf_gamma(r)), lambda(lambda),
                  lambda_pow_r(pow(lambda, r)) { }

        double pdf(double d) const {
                if (d > 0.0) {
                        return 1.0/lambda_pow_r * 1.0/r_gamma * pow(d, r-1.0) * exp(-d/lambda);
                }
                else {
                        return 0.0;
                }
        }
        double rand(const gsl_rng* rng) const {
//                std::cout << gsl_ran_gamma(rng, r, lambda) << std::endl;
                return gsl_ran_gamma(rng, r, lambda);
        }
        double log_gradient(double d) const {
//                return (r-1.0)/d - 1.0/lambda;
                return (1.0/lambda + (r-1.0)/log(1.0-d))/(d-1.0);
        }

        double r;
        double r_gamma;
        double lambda;
        double lambda_pow_r;
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_symbolic_gradient_t : public boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> > {
public:
        pt_symbolic_gradient_t(pt_root_t* root)
                : boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> >() {

                _nodes = root->get_nodes();
                partial_t partial = gradient_rec(root);
                _normalization = poly_sum(partial);

                for (pt_node_t::nodes_t::iterator it = _nodes.begin(); it != _nodes.end(); it++) {
                        const pt_node_t* which = *it;

                        operator[](which) = poly_sum(partial.derivatives[which]);
                }
        }

        polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> normalization() const {
                return _normalization;
        }

        using boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> >::operator[];
        using boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> >::find;

private:
        typedef boost::unordered_map<pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> > derivatives_t;
        typedef boost::array<polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t>, ALPHABET_SIZE+1> carry_t;

        class partial_t : public carry_t {
        public:
                boost::unordered_map<const pt_node_t*, carry_t> derivatives;
        };

        polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> poly_sum(const carry_t& carry) {
                polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> poly_sum;

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> term(1.0);
                        term.exponent()[i] = 1;
                        poly_sum += term*carry[i];
                }
                poly_sum += carry[ALPHABET_SIZE];

                return poly_sum;
        }

        /******************************************************************************
         * Gradient descent
         ******************************************************************************/

        partial_t gradient_leaf(const pt_node_t* node) {
                partial_t partial;
                partial[node->x] += 1.0;
                return partial;
        }
        partial_t gradient_node(
                const pt_node_t* node,
                partial_t& partial_left,
                partial_t& partial_right) {

                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> pm_left (pmut_t(node->left,  true));
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> pm_right(pmut_t(node->right, true));
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> pn_left (pmut_t(node->left,  false));
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> pn_right(pmut_t(node->right, false));

                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> dpm_left ( pn_left );
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> dpm_right( pn_right);
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> dpn_left (-pn_left );
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> dpn_right(-pn_right);

                partial_t partial;
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> poly_sum_left  = poly_sum(partial_left);
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> poly_sum_right = poly_sum(partial_right);

                partial[ALPHABET_SIZE] +=
                        (pn_left *partial_left [ALPHABET_SIZE] + pm_left *poly_sum_left)*
                        (pn_right*partial_right[ALPHABET_SIZE] + pm_right*poly_sum_right);

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        partial[i] += pn_left *pn_right*partial_left [i]*partial_right[ALPHABET_SIZE];
                        partial[i] += pn_left *pm_right*partial_left [i]*poly_sum_right;
                        partial[i] += pn_right*pn_left *partial_right[i]*partial_left [ALPHABET_SIZE];
                        partial[i] += pn_right*pm_left *partial_right[i]*poly_sum_left;
                        partial[i] += pn_left *pn_right*partial_left [i]*partial_right[i];
                }

                for (pt_node_t::nodes_t::iterator it = _nodes.begin(); it != _nodes.end(); it++) {
                        const pt_node_t* which = *it;

                        const polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> deri_sum_left  = poly_sum(partial_left.derivatives[which]);
                        const polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> deri_sum_right = poly_sum(partial_right.derivatives[which]);
                        /* Gradient of sigma
                         */
                        if (node->left == which) {
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        (dpn_left*partial_left [ALPHABET_SIZE]  + dpm_left *poly_sum_left)*
                                        (pn_right*partial_right[ALPHABET_SIZE]  +  pm_right*poly_sum_right);
                        }
                        else if (node->right == which) {
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        ( pn_left *partial_left [ALPHABET_SIZE] +  pm_left *poly_sum_left)*
                                        (dpn_right*partial_right[ALPHABET_SIZE] + dpm_right*poly_sum_right);
                        }
                        else {
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        (pn_left *partial_left [ALPHABET_SIZE]                    + pm_left *poly_sum_left)*
                                        (pn_right*partial_right.derivatives[which][ALPHABET_SIZE] + pm_right*deri_sum_right);
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        (pn_left *partial_left.derivatives [which][ALPHABET_SIZE] + pm_left *deri_sum_left)*
                                        (pn_right*partial_right[ALPHABET_SIZE]                    + pm_right*poly_sum_right);
                        }
                        /* Gradient of phi
                         */
                        if (node->left == which) {
                                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                                        partial.derivatives[which][i] += dpn_left*pn_right*partial_left[ALPHABET_SIZE]*partial_right[i];
                                        partial.derivatives[which][i] += dpm_left*pn_right*poly_sum_left              *partial_right[i];
                                        partial.derivatives[which][i] += dpn_left*pn_right*partial_left[i]            *partial_right[ALPHABET_SIZE];
                                        partial.derivatives[which][i] += dpn_left*pm_right*partial_left[i]            *poly_sum_right;
                                        partial.derivatives[which][i] += dpn_left*pn_right*partial_left[i]            *partial_right[i];
                                }
                        }
                        else if (node->right == which) {
                                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                                        partial.derivatives[which][i] +=  pn_left*dpn_right*partial_left[i]            *partial_right[ALPHABET_SIZE];
                                        partial.derivatives[which][i] +=  pn_left*dpm_right*partial_left[i]            *poly_sum_right;
                                        partial.derivatives[which][i] +=  pn_left*dpn_right*partial_left[ALPHABET_SIZE]*partial_right[i];
                                        partial.derivatives[which][i] +=  pm_left*dpn_right*poly_sum_left              *partial_right[i];
                                        partial.derivatives[which][i] +=  pn_left*dpn_right*partial_left[i]            *partial_right[i];
                                }
                        }
                        else {
                                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                                        partial.derivatives[which][i] += pn_left *pn_right*partial_left[i]*partial_right.derivatives[which][ALPHABET_SIZE];
                                        partial.derivatives[which][i] += pn_left *pm_right*partial_left[i]*deri_sum_right;

                                        partial.derivatives[which][i] += pn_right*pn_left*partial_right[i]*partial_left.derivatives[which] [ALPHABET_SIZE];
                                        partial.derivatives[which][i] += pn_right*pm_left*partial_right[i]*deri_sum_left;

                                        partial.derivatives[which][i] += pn_right*pn_left*partial_right.derivatives[which][i]*partial_left [ALPHABET_SIZE];
                                        partial.derivatives[which][i] += pn_right*pm_left*partial_right.derivatives[which][i]*poly_sum_left;

                                        partial.derivatives[which][i] += pn_left *pn_right*partial_left.derivatives[which][i]*partial_right[ALPHABET_SIZE];
                                        partial.derivatives[which][i] += pn_left *pm_right*partial_left.derivatives[which][i]*poly_sum_right;

                                        partial.derivatives[which][i] += pn_left *pn_right*partial_left [i]*partial_right.derivatives[which][i];
                                        partial.derivatives[which][i] += pn_left *pn_right*partial_right[i]*partial_left.derivatives[which][i];
                                }
                        }
                }
                return partial;
        }

        partial_t gradient_rec(const pt_node_t* node) {
                if (node->leaf()) {
                        return gradient_leaf(node);
                }
                else {
                        partial_t partial_left  = gradient_rec(node->left);
                        partial_t partial_right = gradient_rec(node->right);

                        return gradient_node(node, partial_left, partial_right);
                }
        }

        pt_node_t::nodes_t _nodes;
        polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_tree_t> _normalization;
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_gradient_t : public boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE> > {
public:
        pt_gradient_t(pt_root_t* root)
                : boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE> >() { 

                _nodes = root->get_nodes();
                partial_t partial = gradient_rec(root);
                _normalization = poly_sum(partial);

                for (pt_node_t::nodes_t::iterator it = _nodes.begin(); it != _nodes.end(); it++) {
                        const pt_node_t* which = *it;

                        operator[](which) = poly_sum(partial.derivatives[which]);
                }
        }

        const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& normalization() const {
                return _normalization;
        }

        using boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE> >::operator[];
        using boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE> >::find;

private:
        typedef boost::unordered_map<pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE> > derivatives_t;
        typedef boost::array<polynomial_t<CODE_TYPE, ALPHABET_SIZE>, ALPHABET_SIZE+1> carry_t;

        class partial_t : public carry_t {
        public:
                boost::unordered_map<const pt_node_t*, carry_t> derivatives;
        };

        polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum(const carry_t& carry) {
                polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum;

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> term(1.0);
                        term.exponent()[i] = 1;
                        poly_sum += term*carry[i];
                }
                poly_sum += carry[ALPHABET_SIZE];

                return poly_sum;
        }

        /******************************************************************************
         * Gradient descent
         ******************************************************************************/

        partial_t gradient_leaf(const pt_node_t* node) {
                partial_t partial;
                partial[node->x] += 1.0;
                return partial;
        }
        partial_t gradient_node(
                const pt_node_t* node,
                partial_t& partial_left,
                partial_t& partial_right) {

                double  pm_left  = 1.0-exp(-node->left ->d);
                double  pm_right = 1.0-exp(-node->right->d);
                double  pn_left  = 1.0-pm_left;
                double  pn_right = 1.0-pm_right;

                // double dpm_left  =  pn_left;
                // double dpm_right =  pn_right;
                // double dpn_left  = -pn_left;
                // double dpn_right = -pn_right;

                double dpm_left  =  1.0;
                double dpm_right =  1.0;
                double dpn_left  = -1.0;
                double dpn_right = -1.0;

                partial_t partial;
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum_left  = poly_sum(partial_left);
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum_right = poly_sum(partial_right);

                partial[ALPHABET_SIZE] +=
                        (pn_left *partial_left [ALPHABET_SIZE] + pm_left *poly_sum_left)*
                        (pn_right*partial_right[ALPHABET_SIZE] + pm_right*poly_sum_right);

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        partial[i] += pn_left *pn_right*partial_left [i]*partial_right[ALPHABET_SIZE];
                        partial[i] += pn_left *pm_right*partial_left [i]*poly_sum_right;
                        partial[i] += pn_right*pn_left *partial_right[i]*partial_left [ALPHABET_SIZE];
                        partial[i] += pn_right*pm_left *partial_right[i]*poly_sum_left;
                        partial[i] += pn_left *pn_right*partial_left [i]*partial_right[i];
                }

                for (pt_node_t::nodes_t::iterator it = _nodes.begin(); it != _nodes.end(); it++) {
                        const pt_node_t* which = *it;

                        const polynomial_t<CODE_TYPE, ALPHABET_SIZE> deri_sum_left  = poly_sum(partial_left .derivatives[which]);
                        const polynomial_t<CODE_TYPE, ALPHABET_SIZE> deri_sum_right = poly_sum(partial_right.derivatives[which]);
                        /* Gradient of sigma
                         */
                        if (node->left == which) {
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        (dpn_left*partial_left [ALPHABET_SIZE]  + dpm_left *poly_sum_left)*
                                        (pn_right*partial_right[ALPHABET_SIZE]  +  pm_right*poly_sum_right);
                        }
                        else if (node->right == which) {
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        ( pn_left *partial_left [ALPHABET_SIZE] +  pm_left *poly_sum_left)*
                                        (dpn_right*partial_right[ALPHABET_SIZE] + dpm_right*poly_sum_right);
                        }
                        else {
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        (pn_left *partial_left [ALPHABET_SIZE]                    + pm_left *poly_sum_left)*
                                        (pn_right*partial_right.derivatives[which][ALPHABET_SIZE] + pm_right*deri_sum_right);
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        (pn_left *partial_left.derivatives [which][ALPHABET_SIZE] + pm_left *deri_sum_left)*
                                        (pn_right*partial_right[ALPHABET_SIZE]                    + pm_right*poly_sum_right);
                        }
                        /* Gradient of phi
                         */
                        if (node->left == which) {
                                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                                        partial.derivatives[which][i] += dpn_left*pn_right*partial_left[ALPHABET_SIZE]*partial_right[i];
                                        partial.derivatives[which][i] += dpm_left*pn_right*poly_sum_left              *partial_right[i];
                                        partial.derivatives[which][i] += dpn_left*pn_right*partial_left[i]            *partial_right[ALPHABET_SIZE];
                                        partial.derivatives[which][i] += dpn_left*pm_right*partial_left[i]            *poly_sum_right;
                                        partial.derivatives[which][i] += dpn_left*pn_right*partial_left[i]            *partial_right[i];
                                }
                        }
                        else if (node->right == which) {
                                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                                        partial.derivatives[which][i] +=  pn_left*dpn_right*partial_left[i]            *partial_right[ALPHABET_SIZE];
                                        partial.derivatives[which][i] +=  pn_left*dpm_right*partial_left[i]            *poly_sum_right;
                                        partial.derivatives[which][i] +=  pn_left*dpn_right*partial_left[ALPHABET_SIZE]*partial_right[i];
                                        partial.derivatives[which][i] +=  pm_left*dpn_right*poly_sum_left              *partial_right[i];
                                        partial.derivatives[which][i] +=  pn_left*dpn_right*partial_left[i]            *partial_right[i];
                                }
                        }
                        else {
                                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                                        partial.derivatives[which][i] += pn_left *pn_right*partial_left[i]*partial_right.derivatives[which][ALPHABET_SIZE];
                                        partial.derivatives[which][i] += pn_left *pm_right*partial_left[i]*deri_sum_right;

                                        partial.derivatives[which][i] += pn_right*pn_left*partial_right[i]*partial_left.derivatives[which] [ALPHABET_SIZE];
                                        partial.derivatives[which][i] += pn_right*pm_left*partial_right[i]*deri_sum_left;

                                        partial.derivatives[which][i] += pn_right*pn_left*partial_right.derivatives[which][i]*partial_left [ALPHABET_SIZE];
                                        partial.derivatives[which][i] += pn_right*pm_left*partial_right.derivatives[which][i]*poly_sum_left;

                                        partial.derivatives[which][i] += pn_left *pn_right*partial_left.derivatives[which][i]*partial_right[ALPHABET_SIZE];
                                        partial.derivatives[which][i] += pn_left *pm_right*partial_left.derivatives[which][i]*poly_sum_right;

                                        partial.derivatives[which][i] += pn_left *pn_right*partial_left [i]*partial_right.derivatives[which][i];
                                        partial.derivatives[which][i] += pn_left *pn_right*partial_right[i]*partial_left.derivatives[which][i];
                                }
                        }
                }
                return partial;
        }

        partial_t gradient_rec(const pt_node_t* node) {
                if (node->leaf()) {
                        return gradient_leaf(node);
                }
                else {
                        partial_t partial_left  = gradient_rec(node->left);
                        partial_t partial_right = gradient_rec(node->right);

                        return gradient_node(node, partial_left, partial_right);
                }
        }

        pt_node_t::nodes_t _nodes;
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> _normalization;
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_gradient_ascent_t
{
public:
        pt_gradient_ascent_t(alignment_t<CODE_TYPE>& alignment,
                             const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
                             double r, double lambda, double epsilon = 0.001)
                : alignment(alignment), alpha(alpha), gamma_distribution(r, lambda), epsilon(epsilon) { }

        /* posterior (not normalized) */
        double log_posterior(pt_node_t::nodes_t& nodes) {
                double result = 0;
                // initialize sum with gradients of the gamma distribution
                for (pt_node_t::nodes_t::iterator it = nodes.begin(); it != nodes.end(); it++) {
                        result += log(gamma_distribution.pdf((*it)->d));
                }

                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        const pt_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(alignment.tree);
                        // loop over monomials
                        double tmp = -HUGE_VAL;
                        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator ut = polynomial.begin();
                             ut != polynomial.end(); ut++) {
                                tmp = logadd(tmp, log(ut->coefficient()) + mbeta_log(ut->exponent(), alpha) - mbeta_log(alpha));
                        }
                        result += tmp;
                }
                return result;
        }

        double run() {
                pt_node_t::nodes_t nodes = alignment.tree->get_nodes();
                boost::unordered_map<pt_node_t*, double> sum;
                double total = 0.0;

                // initialize sum with gradients of the gamma distribution
                for (pt_node_t::nodes_t::iterator it = nodes.begin(); it != nodes.end(); it++) {
//                        sum[*it] = gamma_distribution.log_gradient((*it)->d);
                       sum[*it] = gamma_distribution.log_gradient(1.0-exp(-(*it)->d));
                }
                // loop through the alignment
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        pt_gradient_t<CODE_TYPE, ALPHABET_SIZE> gradient(alignment.tree);

                        double norm = 0.0;
                        // loop over monomials
                        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator ut = gradient.normalization().begin();
                             ut != gradient.normalization().end(); ut++) {
                                norm += ut->coefficient()*exp(mbeta_log(ut->exponent(), alpha));
                        }
                        // loop over nodes
                        for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                                double result = 0;
                                // loop over monomials
                                for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator ut = gradient[*is].begin();
                                     ut != gradient[*is].end(); ut++) {
                                        result += ut->coefficient()*exp(mbeta_log(ut->exponent(), alpha));
                                }
                                sum[*is] += result/norm;
                        }
                }
                print_debug("log posterior: %f\n", log_posterior(nodes));
                // apply result
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
//                        double tmp = epsilon/(double)alignment.length*sum[*is];
//                        (*is)->d  += tmp;
                        double tmp = 1.0 - exp(-(*is)->d);
                        tmp       += epsilon/(double)alignment.length*sum[*is];
                        (*is)->d   = -log(1.0-tmp);
                        total     += fabs(tmp);
                }
                return total;
        }
        void run(size_t max, double stop = 0.0) {
                for (size_t i = 0; i < max; i++) {
                        std::stringstream ss;
                        double total = run();
                        alignment.tree->print(ss, true);
                        std::cerr << "total change: " << total << std::endl;
                        std::cerr << ss.str()                  << std::endl;
                        if (total < stop) {
                                break;
                        }
                }
        }

private:
        alignment_t<CODE_TYPE>& alignment;
        exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha;
        gamma_distribution_t gamma_distribution;
        double epsilon;
};

class jumping_distribution_t
{
public:
        virtual double p(double d_old, double d_new) const = 0;
        virtual double sample(double d_old) = 0;
};

class normal_jump_t : public jumping_distribution_t
{
public:
        normal_jump_t(double sigma_square = 0.05)
                : sigma_square(sigma_square) {

                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                _r = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(_r, seed);
        }
        ~normal_jump_t() {
                gsl_rng_free(_r);
        }

        double p(double d_old, double d_new) const {
                return 1.0;
        }
        double sample(double d_old) {
                return old_d+gsl_ran_gaussian(_r, sigma_square);
        }
private:
        double sigma_square;
        gsl_rng * _r;
};

class gamma_jump_t : public jumping_distribution_t
{
public:
        gamma_jump_t(double r, double lambda)
                : gamma_distribution(r, lambda) {

                struct timeval tv;
                gettimeofday(&tv, NULL);
                time_t seed = tv.tv_sec*tv.tv_usec;

                srand(seed);
                _r = gsl_rng_alloc(gsl_rng_default);
                gsl_rng_set(_r, seed);
        }
        ~gamma_jump_t() {
                gsl_rng_free(_r);
        }

        // p(d_new -> d_old)/p(d_old -> d_new) = p(d_old)/p(d_new)
        double p(double d_old, double d_new) const {
                
                return gamma_distribution.pdf(d_old)/gamma_distribution.pdf(d_new);
        }
        double sample(double d_old) {
                return old_d+gsl_ran_gaussian(_r, sigma_square);
        }
private:
        gamma_distribution_t gamma_distribution;
        gsl_rng * _r;
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_metropolis_hastings_t
{
public:
        pt_metropolis_hastings_t(alignment_t<CODE_TYPE>& alignment,
                                 const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
                                 jumping_distribution_t& jumping_distribution)
                : alignment(alignment), alpha(alpha), jumping_distribution(jumping_distribution) {
        }

        double log_likelihood() {
                double result = 0;
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        const pt_polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial(alignment.tree);
                        // loop over monomials
                        double tmp = -HUGE_VAL;
                        for (typename polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator ut = polynomial.begin();
                             ut != polynomial.end(); ut++) {
                                tmp = logadd(tmp, log(ut->coefficient()) + mbeta_log(ut->exponent(), alpha) - mbeta_log(alpha));
                        }
                        result += tmp;
                }
                return result;
        }
        double sample(boost::unordered_map<pt_node_t*, double>& means, size_t i, pt_node_t* node, double log_likelihood_ref) {
                double rho;
                double x;
                double log_likelihood_new;

                double old_d = node->d;
//                double new_d = old_d+gsl_ran_gaussian(_r, sigma_square);
//                double new_d = gamma_distribution.rand(_r);
                double new_d = jumping_distribution.sample(old_d);

                // compute new log likelihood
                node->d            = new_d;
                log_likelihood_new = log_likelihood();

                // compute acceptance probability
//                rho = std::min(1.0, exp(log_likelihood_new-log_likelihood_ref)*gamma_distribution.pdf(new_d)/gamma_distribution.pdf(old_d));
//                rho = std::min(1.0, exp(log_likelihood_new-log_likelihood_ref));
                rho = exp(log_likelihood_new-log_likelihood_ref)*gamma_distribution.pdf(new_d)/gamma_distribution.pdf(old_d)*jumping_distribution.p(d_old, d_new);
                x   = gsl_ran_flat(_r, 0.0, 1.0);
                if (x <= std::min(1.0, rho)) {
                        // sample accepted
                        std::cout << "Sample accepted: " << new_d << std::endl;
                }
                else {
                        // sample rejected
                        std::cout << "Sample rejected: " << new_d << std::endl;
                        node->d = old_d;
                }
                return log_likelihood_new;
        }
        double sample(boost::unordered_map<pt_node_t*, double>& means, size_t i, pt_node_t::nodes_t& nodes) {
                double log_likelihood_ref = log_likelihood();
                // loop over nodes
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                        log_likelihood_ref = sample(means, i, *is, log_likelihood_ref);
                }
                // update means
                std::cout << alignment.tree << std::endl;
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                        means[*is] = ((double)i*means[*is] + (*is)->d)/(double)(i+1);
                        std::cout << "mean: " << means[*is] << std::endl;
                }
                return log_likelihood_ref;
        }
        double sample(size_t burnin, size_t n) {
                boost::unordered_map<pt_node_t*, double> means;
                pt_node_t::nodes_t nodes = alignment.tree->get_nodes();
                double result;
                // burn in
                for (size_t i = 0; i < burnin; i++) {
                        sample(means, i, nodes);
                }
                // sample n times
                for (size_t i = 0; i < n; i++) {
                        result = sample(means, i, nodes);
                }
                // insert means into the tree
                for (boost::unordered_map<pt_node_t*, double>::iterator it = means.begin(); it != means.end(); it++) {
                        it->first->d = it->second;
                }
                return result;
        }

private:
        alignment_t<CODE_TYPE>& alignment;
        exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha;

        jumping_distribution_t& jumping_distribution;
};

#endif /* PHYLOTREE_GRADIENT_HH */
