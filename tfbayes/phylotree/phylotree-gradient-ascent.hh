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

#ifndef PHYLOTREE_GRADIENT_GRADIENT_HH
#define PHYLOTREE_GRADIENT_GRADIENT_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <set>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <algorithm> /* std::min */

#include <tfbayes/phylotree/phylotree-gradient.hh>
#include <tfbayes/phylotree/statistics.hh>

/* This is a gradient ascent method to compute the maximum posterior
 * value with an adaptive step-size similar to resilient
 * backpropagation (Rprop).
 */
template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_gradient_ascent_t
{
public:
        pt_gradient_ascent_t(pt_root_t& tree,
                             alignment_t<CODE_TYPE>& alignment,
                             const exponent_t<CODE_TYPE, ALPHABET_SIZE>& alpha,
                             double r, double lambda,
                             double epsilon = 0.001, double eta = 0.1)
                : tree(tree), alignment(alignment), alpha(alpha), gamma_distribution(r, lambda),
                  epsilon(epsilon), eta(eta) { }

        /* posterior (not normalized) */
        double log_posterior(const pt_node_t::nodes_t& nodes) {
                double result = 0;
                // initialize sum with gradients of the gamma distribution
                for (pt_node_t::nodes_t::const_iterator it = nodes.begin(); it != nodes.end(); it++) {
                        result += log(gamma_distribution.pdf((*it)->d));
                }

                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        const polynomial_t<CODE_TYPE, ALPHABET_SIZE> polynomial = pt_polynomial<CODE_TYPE, ALPHABET_SIZE>(tree, *it);
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

        double run(const pt_node_t::nodes_t& nodes) {
                boost::unordered_map<pt_node_t*, double> sum;
                double total = 0.0;

                // initialize sum with gradients of the gamma distribution
                for (pt_node_t::nodes_t::const_iterator it = nodes.begin(); it != nodes.end(); it++) {
                        sum[*it] = gamma_distribution.log_gradient((*it)->d);
                }
                // loop through the alignment
                for (typename alignment_t<CODE_TYPE>::iterator it = alignment.begin(); it != alignment.end(); it++) {
                        pt_gradient_t<CODE_TYPE, ALPHABET_SIZE> gradient(tree, *it);

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
                        double step;
                        if (sum[*is] > 0) {
                                step =  node_epsilon[*is];
                        }
                        else {
                                step = -node_epsilon[*is];
                        }
                        (*is)->d  = std::max(0.0, (*is)->d+step);
                        total    += fabs(step);
                        if (sum_prev[&**is]*sum[&**is] > 0) {
                                node_epsilon[&**is] *= 1.0+eta;
                        }
                        if (sum_prev[&**is]*sum[&**is] < 0) {
                                node_epsilon[&**is] *= 1.0-eta;
                        }
                }
                sum_prev = sum;

                return total;
        }
        void run(size_t max, double stop = 0.0, bool print = true) {
                pt_node_t::nodes_t nodes = tree.nodes;
                for (pt_node_t::nodes_t::const_iterator is = tree.nodes.begin();
                     is != tree.nodes.end(); is++) {
                        node_epsilon[*is] = epsilon;
                }
                for (size_t i = 0; i < max; i++) {
                        double total = run(nodes);
                        if (print) {
                                std::cout << "total change:  "   << total << std::endl
                                          << newick_format(tree) << std::endl;
                                if (total < stop) {
                                        break;
                                }
                        }
                }
        }

private:
        pt_root_t tree;
        alignment_t<CODE_TYPE> alignment;
        exponent_t<CODE_TYPE, ALPHABET_SIZE> alpha;
        gamma_distribution_t gamma_distribution;
        double epsilon;
        double eta;
        boost::unordered_map<pt_node_t*, double> node_epsilon;
        boost::unordered_map<pt_node_t*, double> sum_prev;
};

#endif /* PHYLOTREE_GRADIENT_GRADIENT_HH */
