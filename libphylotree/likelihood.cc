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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <phylotree.hh>
#include <utility.hh>

using namespace std;

typedef boost::array<polynomial_t<code_t, alphabet_size>, alphabet_size+1> carry_t;

static
polynomial_t<code_t, alphabet_size> poly_sum(const carry_t& carry)
{
        polynomial_t<code_t, alphabet_size> poly_sum;

        for (size_t i = 0; i < alphabet_size; i++) {
                polynomial_term_t<code_t, alphabet_size> term;
                term.exponent()[i] = 1;
                poly_sum += term*carry[i];
        }
        poly_sum += carry[alphabet_size];

        return poly_sum;
}

static
carry_t pt_likelihood_leaf(pt_node_t* node) {
        carry_t carry;
        carry[node->x] += 1;
        return carry;
}

static
carry_t pt_likelihood_node(
        const pt_node_t* node,
        const carry_t& carry_left,
        const carry_t& carry_right)
{
        carry_t carry;
        const polynomial_t<code_t, alphabet_size> poly_sum_left  = poly_sum(carry_left);
        const polynomial_t<code_t, alphabet_size> poly_sum_right = poly_sum(carry_right);

        /*  M, M */
        double p = (1.0-exp(-node->left->d))*(1.0-exp(-node->right->d));
        carry[alphabet_size] += p*poly_sum_left*poly_sum_right;
        for (size_t i = 0; i < alphabet_size+1; i++) {
                /* ¬M, M */
                p = exp(-node->left->d)*(1.0-exp(-node->right->d));
                carry[i] += p*carry_left[i]*poly_sum_right;
                /*  M,¬M */
                p = (1.0-exp(-node->left->d))*exp(-node->right->d);
                carry[i] += p*poly_sum_left*carry_right[i];
                /* ¬M,¬M */
                p = exp(-node->left->d)*exp(-node->right->d);
                carry[i] += p*carry_left[i]*carry_right[i];
                if (i < alphabet_size) {
                        carry[i] += p*carry_left[i]*carry_right[alphabet_size];
                        carry[i] += p*carry_left[alphabet_size]*carry_right[i];
                }
        }

        return carry;
}

static
carry_t pt_likelihood_rec(pt_node_t* node)
{
        if (node->leaf()) {
                return pt_likelihood_leaf(node);
        }
        else {
                const carry_t carry_left  = pt_likelihood_rec(node->left);
                const carry_t carry_right = pt_likelihood_rec(node->right);

                return pt_likelihood_node(node, carry_left, carry_right);
        }
}

polynomial_t<code_t, alphabet_size>
pt_likelihood(pt_root_t* node) {
        return poly_sum(pt_likelihood_rec(node));
}
