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

#ifndef __TFBAYES_PHYLOTREE_PHYLOTREE_GRADIENT_GRADIENT_HH__
#define __TFBAYES_PHYLOTREE_PHYLOTREE_GRADIENT_GRADIENT_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <set>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <algorithm> /* std::min */
#include <limits>

#include <tfbayes/phylotree/phylotree-gradient.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/statistics.hh>

/* This is a gradient ascent method to compute the maximum posterior
 * value with an adaptive step-size similar to resilient
 * backpropagation (Rprop).
 */
template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_gradient_ascent_t
{
public:
        pt_gradient_ascent_t(const pt_root_t& tree,
                             const alignment_map_t<AC>& alignment,
                             const std::vector<double>& alpha,
                             double shape,
                             double scale,
                             thread_pool_t& thread_pool,
                             double epsilon = 0.001,
                             double eta = 0.1)
                : _tree_              (tree),
                  _alignment_         (alignment),
                  _thread_pool_       (thread_pool),
                  _alpha_             (alpha.begin(), alpha.end()),
                  _gamma_distribution_(shape, scale),
                  _epsilon_           (epsilon),
                  _eta_               (eta),
                  node_epsilon        (tree.n_nodes),
                  sum_prev            (tree.n_nodes)
                { }

        double run(const pt_node_t::nodes_t& nodes) {
                std::vector<double> sum(nodes.size());
                double total = 0.0;

                // initialize sum with gradients of the gamma distribution
                for (pt_node_t::nodes_t::const_iterator it = nodes.begin(); it != nodes.end(); it++) {
                        sum[(*it)->id] = _gamma_distribution_.log_gradient((*it)->d);
                }
                // loop through the alignment
                for (typename alignment_map_t<AC>::const_iterator it = _alignment_.begin();
                     it != _alignment_.end(); it++) {
                        // compute the likleihood with gradient information
                        pt_polynomial_derivative_t<AS, PC> likelihood = pt_polynomial_derivative<AS, AC, PC>(_tree_, it->first);

                        double marginal_likelihood = 0.0;
                        // loop over monomials
                        for (typename polynomial_t<AS, PC>::const_iterator ut = likelihood.begin();
                             ut != likelihood.end(); ut++) {
                                marginal_likelihood += ut->coefficient()*exp(mbeta_log(ut->exponent(), _alpha_));
                        }
                        // loop over nodes
                        for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                                const pt_node_t& node = **is;
                                double result = 0;
                                // loop over monomials
                                for (typename polynomial_t<AS, PC>::const_iterator ut = likelihood.derivative()[node.id].begin();
                                     ut != likelihood.derivative()[node.id].end(); ut++) {
                                        result += ut->coefficient()*exp(mbeta_log(ut->exponent(), _alpha_));
                                }
                                sum[node.id] += it->second*result/marginal_likelihood;
                        }
                }
                // apply result
                for (pt_node_t::nodes_t::const_iterator is = nodes.begin(); is != nodes.end(); is++) {
                        pt_node_t& node = **is;
                        double step;
                        if (sum[node.id] > 0) {
                                step =  node_epsilon[node.id];
                        }
                        else {
                                step = -node_epsilon[node.id];
                        }
                        node.d  = std::max(0.0, node.d+step);
                        total  += std::abs(step);
                        if (sum_prev[node.id]*sum[node.id] > 0) {
                                node_epsilon[node.id] *= 1.0+_eta_;
                        }
                        if (sum_prev[node.id]*sum[node.id] < 0) {
                                node_epsilon[node.id] *= 1.0-_eta_;
                        }
                }
                sum_prev = sum;

                return total;
        }
        void run(size_t max, double stop = 0.0, bool verbose = true) {
                pt_node_t::nodes_t nodes = _tree_.nodes;
                for (pt_node_t::nodes_t::const_iterator is = _tree_.nodes.begin();
                     is != _tree_.nodes.end(); is++) {
                        node_epsilon[(*is)->id] = _epsilon_;
                }
                for (size_t i = 0; i < max; i++) {
                        double total = run(nodes);
                        if (verbose) {
                                std::cout << "total change:  "   << total << std::endl
                                          << newick_format(_tree_) << std::endl;
                                if (total < stop) {
                                        break;
                                }
                        }
                }
        }

protected:
        pt_root_t _tree_;
        const alignment_map_t<AC>& _alignment_;
        // a thread pool for computing likelihoods
        thread_pool_t& _thread_pool_;
        exponent_t<AS, PC> _alpha_;
        gamma_distribution_t _gamma_distribution_;
        double _epsilon_;
        double _eta_;
        std::vector<double> node_epsilon;
        std::vector<double> sum_prev;
};

#endif /* __TFBAYES_PHYLOTREE_PHYLOTREE_GRADIENT_GRADIENT_HH__ */
