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

static
void pt_likelihood_leaf(pt_node_t<code_t, alphabet_size>* node) {
        polynomial_term_t<code_t, alphabet_size> term;
        term.exponent()[node->x]   = 1;
        node->poly_sum            += term;
        node->carry[node->x]      += 1;
}

static
void pt_likelihood_node(pt_node_t<code_t, alphabet_size>* node) {
        /*  M, M */
        double p = (1.0-exp(-node->left->d))*(1.0-exp(-node->right->d));
        node->carry[alphabet_size] += p*node->left->poly_sum*node->right->poly_sum;
        for (size_t i = 0; i < alphabet_size+1; i++) {
                /* ¬M, M */
                p = exp(-node->left->d)*(1.0-exp(-node->right->d));
                node->carry[i] += p*node->left->carry[i]*node->right->poly_sum;
                /*  M,¬M */
                p = (1.0-exp(-node->left->d))*exp(-node->right->d);
                node->carry[i] += p*node->left->poly_sum*node->right->carry[i];
                /* ¬M,¬M */
                p = exp(-node->left->d)*exp(-node->right->d);
                node->carry[i] += p*node->left->carry[i]*node->right->carry[i];
                if (i < alphabet_size) {
                        node->carry[i] += p*node->left->carry[i]*node->right->carry[alphabet_size];
                        node->carry[i] += p*node->left->carry[alphabet_size]*node->right->carry[i];
                }
        }

        /* generate sum */
        for (size_t i = 0; i < alphabet_size; i++) {
                polynomial_term_t<code_t, alphabet_size> term;
                term.exponent()[i] = 1;
                node->poly_sum    += term*node->carry[i];
        }
        node->poly_sum += node->carry[alphabet_size];
}

static
void pt_likelihood_rec(pt_node_t<code_t, alphabet_size>* node)
{
        node->init();

        if (node->leaf()) {
                pt_likelihood_leaf(node);
        }
        else {
                pt_likelihood_rec(node->left);
                pt_likelihood_rec(node->right);

                pt_likelihood_node(node);
        }
}

polynomial_t<code_t, alphabet_size>
pt_likelihood(pt_root_t<code_t, alphabet_size>* node) {
        pt_likelihood_rec(node);

        return node->poly_sum;
}
