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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <tfbayes/phylotree/phylotree-simplify.hh>

using namespace std;

incomplete_expression_t pt_simplify_leaf(pt_leaf_t* node) {
        incomplete_expression_t expression;
        incomplete_nodeterm_t term;
        term.incomplete().insert(node);
        expression += term;

        return expression;
}

incomplete_expression_t pt_simplify_root(
        pt_node_t* node,
        const incomplete_expression_t& expression_left,
        const incomplete_expression_t& expression_right) {
        incomplete_expression_t expression1 = expression_left*expression_right;
        incomplete_expression_t expression2;

        for (incomplete_expression_t::const_iterator it = expression1.begin(); it != expression1.end(); it++) {
                incomplete_nodeterm_t term(*it);
                expression2 += term.complete();
        }
        return expression2;
}

incomplete_expression_t pt_simplify_node(
        pt_node_t* node,
        const incomplete_expression_t& expression_left,
        const incomplete_expression_t& expression_right) {
        incomplete_expression_t expression1 = expression_left*expression_right;
        incomplete_expression_t expression2;

        for (incomplete_expression_t::const_iterator it = expression1.begin(); it != expression1.end(); it++) {
                incomplete_nodeterm_t term(*it);
                if (term.incomplete().empty()) {
                        expression2 += term;
                }
                else {
                        expression2 += (1.0-node->mutation_probability())*term;
                        expression2 +=      node->mutation_probability() *term.complete();
                }
        }
        return expression2;
}

incomplete_expression_t pt_simplify_rec(pt_node_t* node) {
        if (node->leaf()) {
                return pt_simplify_leaf(static_cast<pt_leaf_t*>(node));
        }
        else {
                const incomplete_expression_t expression_left  = pt_simplify_rec(node->left);
                const incomplete_expression_t expression_right = pt_simplify_rec(node->right);

                if (node->root()) {
                        return pt_simplify_root(node, expression_left, expression_right);
                }
                else {
                        return pt_simplify_node(node, expression_left, expression_right);
                }
        }
}

incomplete_expression_t pt_simplify(pt_root_t* node) {
        return pt_simplify_rec(node);
}
