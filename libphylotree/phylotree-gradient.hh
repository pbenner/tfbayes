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

#include <set>
#include <boost/unordered_map.hpp>

#include <tfbayes/polynomial.hh>

#include <phylotree.hh>
#include <phylotree-gradient-coefficient.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_gradient_t : public boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> > {
public:
        pt_gradient_t(const pt_root_t* root)
                : boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> >() {

                get_nodes(root);
                partial_t partial = gradient_rec(root);
                _normalization = poly_sum(partial);

                for (nodes_t::iterator it = _nodes.begin(); it != _nodes.end(); it++) {
                        const pt_node_t* which = *it;

                        operator[](which) = poly_sum(partial.derivatives[which]);
                }
        }

        double eval(const pt_node_t* which, const boost::array<double, ALPHABET_SIZE>& p) const {
                return find(which)->second.eval(p)/normalization().eval(p);
        }
        polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> normalization() const {
                return _normalization;
        }

        using boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> >::operator[];
        using boost::unordered_map<const pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> >::find;

private:
        typedef boost::unordered_map<pt_node_t*, polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> > derivatives_t;
        typedef boost::array<polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t>, ALPHABET_SIZE+1> carry_t;

        class partial_t : public carry_t {
        public:
                boost::unordered_map<const pt_node_t*, carry_t> derivatives;
        };

        polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> poly_sum(const carry_t& carry) {
                polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> poly_sum;

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> term(1.0);
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

                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> pm_left (pmut_t(node->left,  true));
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> pm_right(pmut_t(node->right, true));
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> pn_left (pmut_t(node->left,  false));
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> pn_right(pmut_t(node->right, false));

                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> dpm_left ( pn_left );
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> dpm_right( pn_right);
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> dpn_left (-pn_left );
                polynomial_term_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> dpn_right(-pn_right);

                partial_t partial;
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> poly_sum_left  = poly_sum(partial_left);
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> poly_sum_right = poly_sum(partial_right);

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

                for (nodes_t::iterator it = _nodes.begin(); it != _nodes.end(); it++) {
                        const pt_node_t* which = *it;

                        const polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> deri_sum_left  = poly_sum(partial_left.derivatives[which]);
                        const polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> deri_sum_right = poly_sum(partial_right.derivatives[which]);
                        /* Gradient of sigma
                         */
                        if (node->left == which) {
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        (dpn_left*partial_left [ALPHABET_SIZE] + pn_left *poly_sum_left)*
                                        (pn_right*partial_right[ALPHABET_SIZE] + pm_right*poly_sum_right);
                        }
                        else if (node->right == which) {
                                partial.derivatives[which][ALPHABET_SIZE] +=
                                        ( pn_left *partial_left [ALPHABET_SIZE] + pm_left *poly_sum_left)*
                                        (dpn_right*partial_right[ALPHABET_SIZE] + pn_right*poly_sum_right);
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
                                        partial.derivatives[which][i] += dpn_left * pn_right*partial_left [i]*partial_right[ALPHABET_SIZE];
                                        partial.derivatives[which][i] += dpn_left * pm_right*partial_left [i]*poly_sum_right;
                                        partial.derivatives[which][i] +=  pn_right*dpn_left *partial_right[i]*partial_left [ALPHABET_SIZE];
                                        partial.derivatives[which][i] +=  pn_right* pn_left *partial_right[i]*poly_sum_left;
                                        partial.derivatives[which][i] += dpn_left * pn_right*partial_left [i]*partial_right[i];
                                }
                        }
                        else if (node->right == which) {
                                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                                        partial.derivatives[which][i] +=  pn_left *dpn_right*partial_left [i]*partial_right[ALPHABET_SIZE];
                                        partial.derivatives[which][i] +=  pn_left * pn_right*partial_left [i]*poly_sum_right;
                                        partial.derivatives[which][i] += dpn_right* pn_left *partial_right[i]*partial_left [ALPHABET_SIZE];
                                        partial.derivatives[which][i] += dpn_right* pm_left *partial_right[i]*poly_sum_left;
                                        partial.derivatives[which][i] +=  pn_left *dpn_right*partial_left [i]*partial_right[i];
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

        /******************************************************************************
         * Construct a stack of nodes
         ******************************************************************************/

        typedef std::set<const pt_node_t*> nodes_t;

        void get_nodes(const pt_node_t* node) {
                if (!node->root()) {
                        _nodes.insert(node);
                }
                if (!node->leaf()) {
                        get_nodes(node->left);
                        get_nodes(node->right);
                }
        }

        nodes_t _nodes;
        polynomial_t<CODE_TYPE, ALPHABET_SIZE, mutation_coefficient_t> _normalization;
};

#endif /* PHYLOTREE_GRADIENT_HH */
