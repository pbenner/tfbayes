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

ostream& operator<< (ostream& o, const node_set_t& node_set) {
        for (node_set_t::const_iterator it = node_set.begin(); it != node_set.end(); it++) {
                if (it != node_set.begin()) {
                        o << ",";
                }
                o << (*it)->x;
        }

        return o;
}

ostream& operator<< (ostream& o, const incomplete_exponent_t& exponent) {
        for (incomplete_exponent_t::const_iterator it = exponent.begin(); it != exponent.end(); it++) {
                o << "phi("
                  << *it
                  << ")";
        }
        if (!exponent.incomplete().empty()) {
                o << "I("
                  << exponent.incomplete()
                  << ")";
        }
        return o;
}

ostream& operator<< (ostream& o, const incomplete_term_t& term) {
        if (term.coefficient() != 1.0) {
                o << term.coefficient()
                  << " ";
        }
        o << (const incomplete_exponent_t&)term;

        return o;
}

ostream& operator<< (ostream& o, const incomplete_polynomial_t& polynomial) {
        for (incomplete_polynomial_t::const_iterator it = polynomial.begin(); it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + ";
                }
                o << *it;
        }

        return o;
}

size_t hash_value(const node_set_t& set) {
        size_t seed = 0;
        for (node_set_t::const_iterator it = set.begin(); it != set.end(); it++) {
                boost::hash_combine(seed, (void*)*it);
        }
        return seed;
}

size_t hash_value(const incomplete_exponent_t& exponent) {
        size_t seed = 0;
        for (incomplete_exponent_t::const_iterator it = exponent.begin(); it != exponent.end(); it++) {
                boost::hash_combine(seed, hash_value(*it));
        }
        boost::hash_combine(seed, hash_value(exponent.incomplete()));
        return seed;
}

static
double mutation_probability(pt_node_t<code_t, alphabet_size>* node)
{
        return 1.0-exp(-node->d);
}

static
incomplete_polynomial_t
pt_likelihood_leaf(pt_node_t<code_t, alphabet_size>* node) {
        incomplete_polynomial_t poly;
        incomplete_term_t term;
        term.incomplete().insert(node);
        poly += term;

        return poly;
}

static
incomplete_polynomial_t
pt_likelihood_root(
        pt_node_t<code_t, alphabet_size>* node,
        const incomplete_polynomial_t& poly_left,
        const incomplete_polynomial_t& poly_right)
{
        incomplete_polynomial_t poly1 = poly_left*poly_right;
        incomplete_polynomial_t poly2;

        for (incomplete_polynomial_t::const_iterator it = poly1.begin(); it != poly1.end(); it++) {
                incomplete_term_t term(*it);
                poly2 += term.complete();
        }

        return poly2;
}

static
incomplete_polynomial_t
pt_likelihood_node(
        pt_node_t<code_t, alphabet_size>* node,
        const incomplete_polynomial_t& poly_left,
        const incomplete_polynomial_t& poly_right)
{
        incomplete_polynomial_t poly1 = poly_left*poly_right;
        incomplete_polynomial_t poly2;

        for (incomplete_polynomial_t::const_iterator it = poly1.begin(); it != poly1.end(); it++) {
                incomplete_term_t term(*it);

                if (term.incomplete().empty()) {
                        poly2 += term;
                }
                else {
                        poly2 += (1.0-mutation_probability(node))*term;
                        poly2 +=      mutation_probability(node) *term.complete();
                }
        }

        return poly2;
}

static
incomplete_polynomial_t
pt_likelihood_rec(pt_node_t<code_t, alphabet_size>* node)
{
        node->init();

        if (node->leaf()) {
                return pt_likelihood_leaf(node);
        }
        else {
                const incomplete_polynomial_t poly_left  = pt_likelihood_rec(node->left);
                const incomplete_polynomial_t poly_right = pt_likelihood_rec(node->right);

                if (node->root()) {
                        return pt_likelihood_root(node, poly_left, poly_right);
                }
                else {
                        return pt_likelihood_node(node, poly_left, poly_right);
                }
        }
}

polynomial_t<code_t, alphabet_size>
pt_likelihood_prime(pt_root_t<code_t, alphabet_size>* node) {
        cout << pt_likelihood_rec(node);

        return node->poly_sum;
}
