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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/polynomial.hh>

#include <phylotree.hh>
#include <tfbayes/polynomial.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_polynomial_t : public polynomial_t<CODE_TYPE, ALPHABET_SIZE> {
public:
        pt_polynomial_t(const pt_root_t* root)
                : polynomial_t<CODE_TYPE, ALPHABET_SIZE>() {
                likelihood(root);
        }
        pt_polynomial_t(const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly)
                : polynomial_t<CODE_TYPE, ALPHABET_SIZE>(poly) { }

        polynomial_t<CODE_TYPE, ALPHABET_SIZE> likelihood(const pt_root_t* root) {
                operator=(poly_sum(likelihood_rec(root)));
                return *this;
        }

private:
        typedef boost::array<polynomial_t<CODE_TYPE, ALPHABET_SIZE>, ALPHABET_SIZE+1> carry_t;

        polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum(const carry_t& carry) {
                polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum;

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> term;
                        term.exponent()[i] = 1;
                        poly_sum += term*carry[i];
                }
                poly_sum += carry[ALPHABET_SIZE];

                return poly_sum;
        }

        carry_t likelihood_leaf(const pt_node_t* node) {
                carry_t carry;
                carry[node->x] += 1;
                return carry;
        }
        carry_t likelihood_node(
                const pt_node_t* node,
                const carry_t& carry_left,
                const carry_t& carry_right) {
                double p;
                carry_t carry;
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum_left  = poly_sum(carry_left);
                const polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_sum_right = poly_sum(carry_right);

                /*  M, M */
                p = (1.0-exp(-node->left->d))*(1.0-exp(-node->right->d));
                carry[ALPHABET_SIZE] += p*poly_sum_left*poly_sum_right;
                /* ¬M, M */
                p = exp(-node->left->d)*(1.0-exp(-node->right->d));
                carry[ALPHABET_SIZE] += p*carry_left[ALPHABET_SIZE]*poly_sum_right;
                /*  M,¬M */
                p = (1.0-exp(-node->left->d))*exp(-node->right->d);
                carry[ALPHABET_SIZE] += p*poly_sum_left*carry_right[ALPHABET_SIZE];
                /* ¬M,¬M */
                p = exp(-node->left->d)*exp(-node->right->d);
                carry[ALPHABET_SIZE] += p*carry_left[ALPHABET_SIZE]*carry_right[ALPHABET_SIZE];

                for (size_t i = 0; i < ALPHABET_SIZE; i++) {
                        /*  M, M */
                        carry[i] += p*carry_left[i]*carry_right[ALPHABET_SIZE];
                        carry[i] += p*carry_left[ALPHABET_SIZE]*carry_right[i];
                        /* ¬M, M */
                        p = exp(-node->left->d)*(1.0-exp(-node->right->d));
                        carry[i] += p*carry_left[i]*poly_sum_right;
                        /*  M,¬M */
                        p = (1.0-exp(-node->left->d))*exp(-node->right->d);
                        carry[i] += p*poly_sum_left*carry_right[i];
                        /* ¬M,¬M */
                        p = exp(-node->left->d)*exp(-node->right->d);
                        carry[i] += p*carry_left[i]*carry_right[i];
                }
                return carry;
        }

        carry_t likelihood_rec(const pt_node_t* node) {
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

#endif /* PHYLOTREE_POLYNOMIAL_HH */
