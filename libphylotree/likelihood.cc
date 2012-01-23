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
        term[node->x] = 1;
        node->poly_sum            += term;
        node->carry[node->x]      += 1;
        node->applicable[node->x]  = true;
}

static
void pt_likelihood_node(pt_node_t<code_t, alphabet_size>* node) {
        /*  M, M */
        const double p = (1.0-exp(-node->left->d))*(1.0-exp(-node->right->d));
        node->carry[alphabet_size] += p*node->left->poly_sum*node->right->poly_sum;
        for (size_t i = 0; i < alphabet_size+1; i++) {
                /* ¬M, M */
                if (node->left->applicable[i]) {
                        const double p = exp(-node->left->d)*(1.0-exp(-node->right->d));
                        node->carry[i] += p*node->left->carry[i]*node->right->poly_sum;
                        node->applicable[i] = true;
                }
                /*  M,¬M */
                if (node->right->applicable[i]) {
                        const double p = (1.0-exp(-node->left->d))*exp(-node->right->d);
                        node->carry[i] += p*node->left->poly_sum*node->right->carry[i];
                        node->applicable[i] = true;
                }
                /* ¬M,¬M */
                if (node->left->applicable[i] && node->right->applicable[i]) {
                        const double p = exp(-node->left->d)*exp(-node->right->d);
                        node->carry[i] += p*node->left->carry[i]*node->right->carry[i];
                        node->applicable[i] = true;
                }
        }

        /* generate sum */
        for (size_t i = 0; i < alphabet_size+1; i++) {
                /* if this is the root then multiply by Px */
                if (node->root() && i < alphabet_size) {
                        polynomial_term_t<code_t, alphabet_size> term;
                        term[i] = 1;
                        node->poly_sum += term*node->carry[i];
                }
                /* otherwise just sum up */
                else {
                        node->poly_sum += node->carry[i];
                }
        }
}

static
void pt_likelihood_rec(pt_node_t<code_t, alphabet_size>* node) {

        if (node->leaf()) {
                pt_likelihood_leaf(node);
        }
        else {
                pt_likelihood_rec(node->left);
                pt_likelihood_rec(node->right);

                pt_likelihood_node(node);
        }
}

double pt_likelihood() {

        pt_leaf_t<code_t, alphabet_size> n2(1, 1.0);
        pt_leaf_t<code_t, alphabet_size> n3(2, 2.0);
        pt_root_t<code_t, alphabet_size> n1(-1, &n2, &n3);

        pt_likelihood_rec(&n1);
        cout << n1.poly_sum << endl;

        boost::array<double, alphabet_size> val;
        val[0] = 0.1;
        val[1] = 0.3;
        val[2] = 0.5;
        val[3] = 0.1;
        cout << "eval: "
             << n1.poly_sum.eval(val)
             << endl;

        return 0.0;
}
