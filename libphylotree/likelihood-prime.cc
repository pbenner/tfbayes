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
#include <vector>

#include <incomplete-polynomial.hh>
#include <phylotree.hh>
#include <utility.hh>

using namespace std;

size_t hash_value(const node_set_t set) {
        size_t seed = 0;
        for (node_set_t::const_iterator it = set.begin(); it != set.end(); it++) {
                boost::hash_combine(seed, (void*)*it);
        }
        return seed;
}

size_t hash_value(const incomplete_term_t term) {
        size_t seed = 0;
        boost::hash_combine(seed, hash_value(term.incomplete));
        return seed;
}

static
void pt_likelihood_leaf(pt_node_t<code_t, alphabet_size>* node) {
        incomplete_polynomial_t poly;
        incomplete_term_t term;
        term.incomplete.push_back(node);
        poly[term] = 1;
}

static
void pt_likelihood_node(pt_node_t<code_t, alphabet_size>* node) {
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
pt_likelihood_prime(pt_root_t<code_t, alphabet_size>* node) {
        pt_likelihood_rec(node);

        return node->poly_sum;
}
