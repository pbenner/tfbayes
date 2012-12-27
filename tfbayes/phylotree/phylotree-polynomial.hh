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

#ifndef PHYLOTREE_POLYNOMIAL_HH
#define PHYLOTREE_POLYNOMIAL_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/utility/polynomial.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_polynomial_class {
public:
        static polynomial_t<CODE_TYPE, ALPHABET_SIZE> likelihood(const pt_root_t* root) {
                return poly_sum(likelihood_rec(root));
        }

private:
        /******************************************************************************
         * typedefs
         ******************************************************************************/

        typedef boost::array<polynomial_t<CODE_TYPE, ALPHABET_SIZE>, ALPHABET_SIZE+1> carry_t;

        /******************************************************************************
         * Algorithm by PY and PB
         ******************************************************************************/

        static polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum(const carry_t& carry) {
                polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum;

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> term(1.0);
                        term.exponent()[i] = 1;
                        poly_sum += term*carry[i];
                }
                poly_sum += carry[ALPHABET_SIZE];

                return poly_sum;
        }

        static carry_t likelihood_leaf(const pt_node_t* node) {
                carry_t carry;
                if (node->x != -1) {
                        // leaf is present and has a nucleotide
                        carry[node->x] += 1.0;
                }
                else {
                        // this leaf should be ignored, no nucleotide
                        // is present
                        carry[ALPHABET_SIZE] = 1.0;
                }
                return carry;
        }
        static carry_t likelihood_node(
                const pt_node_t* node,
                const carry_t& carry_left,
                const carry_t& carry_right) {
                double pm_left  = 1.0-exp(-node->left->d);
                double pm_right = 1.0-exp(-node->right->d);
                carry_t carry;
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum_left  = poly_sum(carry_left);
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum_right = poly_sum(carry_right);

                carry[ALPHABET_SIZE] +=
                        ((1.0-pm_left )*carry_left [ALPHABET_SIZE] + pm_left *poly_sum_left)*
                        ((1.0-pm_right)*carry_right[ALPHABET_SIZE] + pm_right*poly_sum_right);

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        carry[i] += (1.0-pm_left )*(1.0-pm_right)*carry_left [i]*carry_right[ALPHABET_SIZE];
                        carry[i] += (1.0-pm_left )*     pm_right *carry_left [i]*poly_sum_right;
                        carry[i] += (1.0-pm_right)*(1.0-pm_left )*carry_right[i]*carry_left [ALPHABET_SIZE];
                        carry[i] += (1.0-pm_right)*     pm_left  *carry_right[i]*poly_sum_left;
                        carry[i] += (1.0-pm_left )*(1.0-pm_right)*carry_left [i]*carry_right[i];
                }
                return carry;
        }

        static carry_t likelihood_rec(const pt_node_t* node) {
                if (node->leaf()) {
                        return likelihood_leaf(node);
                }
                else {
                        const carry_t carry_left  = likelihood_rec(node->left);
                        const carry_t carry_right = likelihood_rec(node->right);

                        return likelihood_node(node, carry_left, carry_right);
                }
        }
};

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_polynomial(const pt_root_t* root)
{
        return pt_polynomial_class<CODE_TYPE, ALPHABET_SIZE>::likelihood(root);
}

#endif /* PHYLOTREE_POLYNOMIAL_HH */