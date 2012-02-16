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
                seed += (size_t)*it;
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
double mutation_probability(pt_node_t* node)
{
        return 1.0-exp(-node->d);
}
static
double mutation_probability(double d)
{
        return 1.0-exp(-d);
}

static
incomplete_polynomial_t
pt_likelihood_leaf(pt_node_t* node) {
        incomplete_polynomial_t poly;
        incomplete_term_t term;
        term.incomplete().insert(node);
        poly += term;

        return poly;
}

static
incomplete_polynomial_t
pt_likelihood_root(
        pt_node_t* node,
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
        pt_node_t* node,
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
pt_likelihood_rec(pt_node_t* node)
{
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

static
polynomial_term_t<code_t, alphabet_size>
nucleotide_probability(code_t x) {
        polynomial_term_t<code_t, alphabet_size> px;
        px.exponent()[x] = 1;
        return px;
}

static
polynomial_t<code_t, alphabet_size>
mutation_model(code_t x, code_t y, double d) {
        polynomial_t<code_t, alphabet_size> poly;
        polynomial_term_t<code_t, alphabet_size> px = nucleotide_probability(x);
        poly += mutation_probability(d)*px;
        if (x == y) {
                poly += (1-mutation_probability(d));
        }
        return poly;
}

polynomial_t<code_t, alphabet_size>
expand(const node_set_t& node_set) {
        vector<bool> applicable(alphabet_size, false);
        polynomial_t<code_t, alphabet_size> result;
        polynomial_t<code_t, alphabet_size> remainder(1.0);

        for (node_set_t::const_iterator it = node_set.begin(); it != node_set.end(); it++) {
                applicable[(*it)->x] = true;
        }
        for (code_t y = 0; y < alphabet_size; y++) {
                if (applicable[y]) {
                        polynomial_term_t<code_t, alphabet_size> px = nucleotide_probability(y);
                        polynomial_t<code_t, alphabet_size> tmp(px);
                        remainder -= px;
                        for (node_set_t::const_iterator it = node_set.begin(); it != node_set.end(); it++) {
                                tmp *= mutation_model((*it)->x, y, (*it)->d);
                        }
                        result += tmp;
                }
        }
        if (remainder.size() < alphabet_size) {
                polynomial_t<code_t, alphabet_size> tmp(remainder);
                for (node_set_t::const_iterator it = node_set.begin(); it != node_set.end(); it++) {
                        tmp *= mutation_model((*it)->x, alphabet_size, (*it)->d);
                }
                result += tmp;
        }
        return result;
}

polynomial_t<code_t, alphabet_size>
expand(const incomplete_term_t& term) {
        polynomial_t<code_t, alphabet_size> result(1.0);
        for (incomplete_term_t::const_iterator it = term.begin(); it != term.end(); it++) {
                result *= expand(*it);
        }
        return term.coefficient()*result;
}

polynomial_t<code_t, alphabet_size>
expand(const incomplete_polynomial_t& poly) {
        polynomial_t<code_t, alphabet_size> result;
        for (incomplete_polynomial_t::const_iterator it = poly.begin(); it != poly.end(); it++) {
                result += expand(*it);
        }

        return result;
}

polynomial_t<code_t, alphabet_size>
pt_likelihood_prime(pt_root_t* node) {
        polynomial_t<code_t, alphabet_size> poly_sum;
        incomplete_polynomial_t incomplete_poly = pt_likelihood_rec(node);

        return expand(incomplete_poly);
}
