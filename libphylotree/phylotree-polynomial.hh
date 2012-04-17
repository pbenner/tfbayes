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
#include <incomplete-polynomial.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
class pt_polynomial_t : public polynomial_t<CODE_TYPE, ALPHABET_SIZE> {
public:
        pt_polynomial_t(const pt_root_t* root)
                : polynomial_t<CODE_TYPE, ALPHABET_SIZE>() {
                likelihood_py(root);
        }
        pt_polynomial_t(const polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly)
                : polynomial_t<CODE_TYPE, ALPHABET_SIZE>(poly) { }

        polynomial_t<CODE_TYPE, ALPHABET_SIZE> likelihood_py(const pt_root_t* root) {
                operator=(poly_sum(likelihood_py_rec(root)));
                return *this;
        }

        polynomial_t<CODE_TYPE, ALPHABET_SIZE> likelihood(const pt_root_t* root) {
                operator=(likelihood_rec(root).complete());
                return *this;
        }

private:

        /******************************************************************************
         * Algorithm by PY and ST
         ******************************************************************************/
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

        carry_t likelihood_py_leaf(const pt_node_t* node) {
                carry_t carry;
                carry[node->x] += 1;
                return carry;
        }
        carry_t likelihood_py_node(
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

        carry_t likelihood_py_rec(const pt_node_t* node) {
                if (node->leaf()) {
                        return likelihood_py_leaf(node);
                }
                else {
                        const carry_t carry_left  = likelihood_py_rec(node->left);
                        const carry_t carry_right = likelihood_py_rec(node->right);

                        return likelihood_py_node(node, carry_left, carry_right);
                }
        }

        /******************************************************************************
         * Algorithm by PB
         ******************************************************************************/

        incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> likelihood_leaf(const pt_node_t* node) {
                incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly;
                incomplete_term_t<CODE_TYPE, ALPHABET_SIZE> term;
                term.incomplete().insert(node);
                poly += term;

                return poly;
        }

        incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> likelihood_root(
                const pt_node_t* node,
                const incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly_left,
                const incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly_right) {
                incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly1 = poly_left*poly_right;
                incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly2;

                for (typename incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = poly1.begin(); it != poly1.end(); it++) {
                        incomplete_term_t<CODE_TYPE, ALPHABET_SIZE> term(*it);
                        poly2 += term.complete();
                }
                return poly2;
        }

        incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> likelihood_node(
                const pt_node_t* node,
                const incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly_left,
                const incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE>& poly_right) {
                incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly1 = poly_left*poly_right;
                incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly2;

                for (typename incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE>::const_iterator it = poly1.begin(); it != poly1.end(); it++) {
                        incomplete_term_t<CODE_TYPE, ALPHABET_SIZE> term(*it);
                        if (term.incomplete().empty()) {
                                poly2 += term;
                        }
                        else {
                                poly2 += (1.0-node->mutation_probability())*term;
                                poly2 +=      node->mutation_probability() *term.complete();
                        }
                }
                std::cerr << "poly2: " << poly2 << std::endl;
                return poly2;
        }

        incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> likelihood_rec(const pt_node_t* node) {
                if (node->leaf()) {
                        return likelihood_leaf(node);
                }
                else {
                        const incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_left  = likelihood_rec(node->left);
                        const incomplete_polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly_right = likelihood_rec(node->right);

                        if (node->root()) {
                                return likelihood_root(node, poly_left, poly_right);
                        }
                        else {
                                return likelihood_node(node, poly_left, poly_right);
                        }
                }
        }
};

#endif /* PHYLOTREE_POLYNOMIAL_HH */
