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

#ifndef __TFBAYES_PHYLOTREE_PHYLOTREE_GRADIENT_HH__
#define __TFBAYES_PHYLOTREE_PHYLOTREE_GRADIENT_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <set>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <algorithm> /* std::min */

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/distribution.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/polynomial.hh>

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_gradient_t : public boost::unordered_map<pt_node_t::id_t, polynomial_t<AS, PC> > {
public:
        typedef boost::unordered_map<pt_node_t::id_t, polynomial_t<AS, PC> > base_t;

        pt_gradient_t(const pt_root_t& root, const std::vector<AC>& observations)
                : base_t() { 

                _nodes = root.nodes;
                partial_t partial = gradient_rec(root, observations);
                _likelihood = poly_sum(partial);

                for (pt_node_t::nodes_t::const_iterator it = _nodes.begin(); it != _nodes.end(); it++) {
                        const pt_node_t::id_t id = (*it)->id;

                        base_t::operator[](id) = poly_sum(partial.derivatives[id]);
                }
        }

        const polynomial_t<AS, PC>& likelihood() const {
                return _likelihood;
        }

protected:
        typedef boost::unordered_map<pt_node_t::id_t, polynomial_t<AS, PC> > derivatives_t;
        typedef boost::array<polynomial_t<AS, PC>, AS+1> carry_t;

        class partial_t : public carry_t {
        public:
                boost::unordered_map<pt_node_t::id_t, carry_t> derivatives;
        };

        polynomial_t<AS, PC> poly_sum(const carry_t& carry) {
                polynomial_t<AS, PC> poly_sum;

                for (size_t i = 0; i < AS; i++) {
                        polynomial_term_t<AS, PC> term(1.0);
                        term.exponent()[i] = 1;
                        poly_sum += term*carry[i];
                }
                poly_sum += carry[AS];

                return poly_sum;
        }

        /******************************************************************************
         * Gradient descent
         ******************************************************************************/

        partial_t gradient_leaf(
                const pt_leaf_t& leaf,
                const std::vector<AC>& observations) {
                partial_t partial;
                partial[observations[leaf.id]] += 1.0;
                return partial;
        }
        partial_t gradient_node(
                const pt_node_t& node,
                partial_t& partial_left,
                partial_t& partial_right,
                const std::vector<AC>& observations) {

                double  pm_left  = 1.0-std::exp(-node.left ().d);
                double  pm_right = 1.0-std::exp(-node.right().d);
                double  pn_left  = 1.0-pm_left;
                double  pn_right = 1.0-pm_right;

                double dpm_left  =  pn_left;
                double dpm_right =  pn_right;
                double dpn_left  = -pn_left;
                double dpn_right = -pn_right;

                partial_t partial;
                const polynomial_t<AS, PC> poly_sum_left  = poly_sum(partial_left);
                const polynomial_t<AS, PC> poly_sum_right = poly_sum(partial_right);

                partial[AS] +=
                        (pn_left *partial_left [AS] + pm_left *poly_sum_left)*
                        (pn_right*partial_right[AS] + pm_right*poly_sum_right);

                for (size_t i = 0; i < AS; i++) {
                        partial[i] += pn_left *pn_right*partial_left [i]*partial_right[AS];
                        partial[i] += pn_left *pm_right*partial_left [i]*poly_sum_right;
                        partial[i] += pn_right*pn_left *partial_right[i]*partial_left [AS];
                        partial[i] += pn_right*pm_left *partial_right[i]*poly_sum_left;
                        partial[i] += pn_left *pn_right*partial_left [i]*partial_right[i];
                }

                for (pt_node_t::nodes_t::iterator it = _nodes.begin(); it != _nodes.end(); it++) {
                        const pt_node_t::id_t id = (*it)->id;

                        const polynomial_t<AS, PC> deri_sum_left  = poly_sum(partial_left .derivatives[id]);
                        const polynomial_t<AS, PC> deri_sum_right = poly_sum(partial_right.derivatives[id]);
                        /* Gradient of sigma
                         */
                        if (node.left().id == id) {
                                partial.derivatives[id][AS] +=
                                        (dpn_left*partial_left [AS]  + dpm_left *poly_sum_left)*
                                        (pn_right*partial_right[AS]  +  pm_right*poly_sum_right);
                        }
                        else if (node.right().id == id) {
                                partial.derivatives[id][AS] +=
                                        ( pn_left *partial_left [AS] +  pm_left *poly_sum_left)*
                                        (dpn_right*partial_right[AS] + dpm_right*poly_sum_right);
                        }
                        else {
                                partial.derivatives[id][AS] +=
                                        (pn_left *partial_left [AS]                    + pm_left *poly_sum_left)*
                                        (pn_right*partial_right.derivatives[id][AS] + pm_right*deri_sum_right);
                                partial.derivatives[id][AS] +=
                                        (pn_left *partial_left.derivatives [id][AS] + pm_left *deri_sum_left)*
                                        (pn_right*partial_right[AS]                    + pm_right*poly_sum_right);
                        }
                        /* Gradient of phi
                         */
                        if (node.left().id == id) {
                                for (size_t i = 0; i < AS; i++) {
                                        partial.derivatives[id][i] += dpn_left*pn_right*partial_left[AS]*partial_right[ i];
                                        partial.derivatives[id][i] += dpm_left*pn_right*poly_sum_left   *partial_right[ i];
                                        partial.derivatives[id][i] += dpn_left*pn_right*partial_left[ i]*partial_right[AS];
                                        partial.derivatives[id][i] += dpn_left*pm_right*partial_left[ i]*poly_sum_right;
                                        partial.derivatives[id][i] += dpn_left*pn_right*partial_left[ i]*partial_right[ i];
                                }
                        }
                        else if (node.right().id == id) {
                                for (size_t i = 0; i < AS; i++) {
                                        partial.derivatives[id][i] +=  pn_left*dpn_right*partial_left[ i]*partial_right[AS];
                                        partial.derivatives[id][i] +=  pn_left*dpm_right*partial_left[ i]*poly_sum_right;
                                        partial.derivatives[id][i] +=  pn_left*dpn_right*partial_left[AS]*partial_right[ i];
                                        partial.derivatives[id][i] +=  pm_left*dpn_right*poly_sum_left   *partial_right[ i];
                                        partial.derivatives[id][i] +=  pn_left*dpn_right*partial_left[ i]*partial_right[ i];
                                }
                        }
                        else {
                                for (size_t i = 0; i < AS; i++) {
                                        partial.derivatives[id][i] += pn_left *pn_right*partial_left[i]*partial_right.derivatives[id][AS];
                                        partial.derivatives[id][i] += pn_left *pm_right*partial_left[i]*deri_sum_right;

                                        partial.derivatives[id][i] += pn_right*pn_left*partial_right[i]*partial_left .derivatives[id][AS];
                                        partial.derivatives[id][i] += pn_right*pm_left*partial_right[i]*deri_sum_left;

                                        partial.derivatives[id][i] += pn_right*pn_left*partial_right.derivatives[id][i]*partial_left[AS];
                                        partial.derivatives[id][i] += pn_right*pm_left*partial_right.derivatives[id][i]*poly_sum_left;

                                        partial.derivatives[id][i] += pn_left *pn_right*partial_left.derivatives[id][i]*partial_right[AS];
                                        partial.derivatives[id][i] += pn_left *pm_right*partial_left.derivatives[id][i]*poly_sum_right;

                                        partial.derivatives[id][i] += pn_left *pn_right*partial_left [i]*partial_right.derivatives[id][i];
                                        partial.derivatives[id][i] += pn_left *pn_right*partial_right[i]*partial_left .derivatives[id][i];
                                }
                        }
                }
                return partial;
        }

        partial_t gradient_rec(
                const pt_node_t& node,
                const std::vector<AC>& observations) {
                if (node.leaf()) {
                        return gradient_leaf(static_cast<const pt_leaf_t&>(node), observations);
                }
                else {
                        partial_t partial_left  = gradient_rec(node.left (), observations);
                        partial_t partial_right = gradient_rec(node.right(), observations);

                        return gradient_node(node, partial_left, partial_right, observations);
                }
        }

        pt_node_t::nodes_t _nodes;
        polynomial_t<AS, PC> _likelihood;
};

#endif /* __TFBAYES_PHYLOTREE_PHYLOTREE_GRADIENT_HH__ */
