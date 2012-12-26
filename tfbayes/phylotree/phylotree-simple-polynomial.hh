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

#ifndef PHYLOTREE_SIMPLE_POLYNOMIAL_HH
#define PHYLOTREE_SIMPLE_POLYNOMIAL_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-expand.hh>
#include <tfbayes/utility/polynomial.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_simple_polynomial_t : public polynomial_t<CODE_TYPE, ALPHABET_SIZE> {
public:
        pt_simple_polynomial_t(pt_root_t* root)
                : polynomial_t<CODE_TYPE, ALPHABET_SIZE>() {
                likelihood_py(root);
        }

        pt_simple_polynomial_t(const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly)
                : polynomial_t<CODE_TYPE, ALPHABET_SIZE>(poly) { }

        polynomial_t<CODE_TYPE, ALPHABET_SIZE>& likelihood_py(pt_root_t* root) {
                this->operator=(likelihood_rec(root)[ALPHABET_SIZE]);
                return *this;
        }

private:
        typedef boost::array<polynomial_t<CODE_TYPE, ALPHABET_SIZE>, ALPHABET_SIZE+1> carry_t;

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

        carry_t likelihood_leaf(pt_node_t* node) {
                carry_t carry;
                if (node->root()) {
                        carry[ALPHABET_SIZE] = nucleotide_probability<CODE_TYPE, ALPHABET_SIZE>(node->x);
                }
                else {
                        // i: nucleotide of ancestor
                        for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                                carry[i] = mutation_model<CODE_TYPE, ALPHABET_SIZE>(node, i);
                        }
                }
                return carry;
        }
        carry_t likelihood_node(
                pt_node_t* node,
                const carry_t& carry_left,
                const carry_t& carry_right) {
                carry_t carry;

                // i: nucleotide of ancestor
                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        // j: nucleotide of this node
                        for (size_t j = 0; j < ALPHABET_SIZE; j++) {
                                node->x = j;
                                carry[i] += carry_left[j]*carry_right[j]*mutation_model<CODE_TYPE, ALPHABET_SIZE>(node, i);
                        }
                        node->x = -1;
                }
                return carry;
        }
        carry_t likelihood_root(
                pt_node_t* node,
                const carry_t& carry_left,
                const carry_t& carry_right) {
                carry_t carry;

                // j: nucleotide of this node
                for (size_t j = 0; j < ALPHABET_SIZE; j++) {
                        carry[ALPHABET_SIZE] += carry_left[j]*carry_right[j]*nucleotide_probability<CODE_TYPE, ALPHABET_SIZE>(j);
                }
                return carry;
        }

        carry_t likelihood_rec(pt_node_t* node) {
                if (node->leaf()) {
                        return likelihood_leaf(node);
                }
                else {
                        const carry_t carry_left  = likelihood_rec(node->left);
                        const carry_t carry_right = likelihood_rec(node->right);

                        if (node->root()) {
                                return likelihood_root(node, carry_left, carry_right);
                        }
                        else {
                                return likelihood_node(node, carry_left, carry_right);
                        }
                }
        }
};

#endif /* PHYLOTREE_SIMPLE_POLYNOMIAL_HH */
