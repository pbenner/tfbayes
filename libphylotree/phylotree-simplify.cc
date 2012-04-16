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

#include <phylotree-simplify.hh>

using namespace std;

incomplete_polynomial_t pt_simplify_leaf(pt_node_t* node) {
        incomplete_polynomial_t poly;
        incomplete_term_t term;
        term.incomplete().insert(node);
        poly += term;

        return poly;
}

incomplete_polynomial_t pt_simplify_root(
        pt_node_t* node,
        const incomplete_polynomial_t& poly_left,
        const incomplete_polynomial_t& poly_right) {
        incomplete_polynomial_t poly1 = poly_left*poly_right;
        incomplete_polynomial_t poly2;

        for (incomplete_polynomial_t::const_iterator it = poly1.begin(); it != poly1.end(); it++) {
                incomplete_term_t term(*it);
                poly2 += term.complete();
        }
        return poly2;
}

incomplete_polynomial_t pt_simplify_node(
        pt_node_t* node,
        const incomplete_polynomial_t& poly_left,
        const incomplete_polynomial_t& poly_right) {
        incomplete_polynomial_t poly1 = poly_left*poly_right;
        incomplete_polynomial_t poly2;

        for (incomplete_polynomial_t::const_iterator it = poly1.begin(); it != poly1.end(); it++) {
                incomplete_term_t term(*it);
                if (term.incomplete().empty()) {
                        poly2 += term;
                }
                else {
                        poly2 += (1.0-node->mutation_probability())*term;
                        poly2 +=      node->mutation_probability() *term.complete();
                }
        }
        return poly2;
}

incomplete_polynomial_t pt_simplify_rec(pt_node_t* node) {
        if (node->leaf()) {
                return pt_simplify_leaf(node);
        }
        else {
                const incomplete_polynomial_t poly_left  = pt_simplify_rec(node->left);
                const incomplete_polynomial_t poly_right = pt_simplify_rec(node->right);

                if (node->root()) {
                        return pt_simplify_root(node, poly_left, poly_right);
                }
                else {
                        return pt_simplify_node(node, poly_left, poly_right);
                }
        }
}

incomplete_polynomial_t pt_simplify(pt_root_t* node) {
        return pt_simplify_rec(node);
}
