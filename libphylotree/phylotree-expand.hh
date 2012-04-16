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

#ifndef PHYLOTREE_EXPAND_HH
#define PHYLOTREE_EXPAND_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <tfbayes/polynomial.hh>

#include <phylotree.hh>
#include <incomplete-polynomial.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> nucleotide_probability(CODE_TYPE x) {
        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> px(1.0);
        px.exponent()[x] = 1;
        return px;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> mutation_model(const pt_node_t* u, CODE_TYPE y) {
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> poly;
        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> px = nucleotide_probability<CODE_TYPE, ALPHABET_SIZE>(u->x);
        poly += u->mutation_probability()*px;
        if (u->x == y) {
                poly += (1.0-u->mutation_probability());
        }
        return poly;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_expand_rec(
        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> term,
        node_set_t::const_iterator it,
        node_set_t::const_iterator end,
        CODE_TYPE condition)
{
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> result(0.0);
        if(it == end) {
                /* leaf */
                if (condition != ALPHABET_SIZE) {
                        term *= nucleotide_probability<CODE_TYPE, ALPHABET_SIZE>(condition);
                }
                result += term;
        }
        else {
                const double pm   = (*it)->mutation_probability();
                const CODE_TYPE x = (*it)->x; it++;

                /* no mutation */
                if (condition == ALPHABET_SIZE || condition == x) {
                        result += pt_expand_rec((1.0-pm)*term, it, end, x);
                }
                /* mutation */
                term   *= nucleotide_probability<CODE_TYPE, ALPHABET_SIZE>(x);
                result += pt_expand_rec(pm*term, it, end, condition);
        }
        return result;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_expand(const node_set_t& node_set) {
        return pt_expand_rec<CODE_TYPE, ALPHABET_SIZE>(1.0, node_set.begin(), node_set.end(), ALPHABET_SIZE);
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_expand(const incomplete_term_t& term) {
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> result(1.0);
        for (incomplete_term_t::const_iterator it = term.begin(); it != term.end(); it++) {
                result *= pt_expand<CODE_TYPE, ALPHABET_SIZE>(*it);
        }
        return term.coefficient()*result;
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_expand(const incomplete_polynomial_t& poly) {
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> result;
        for (incomplete_polynomial_t::const_iterator it = poly.begin(); it != poly.end(); it++) {
                result += pt_expand<CODE_TYPE, ALPHABET_SIZE>(*it);
        }

        return result;
}

#endif /* PHYLOTREE_EXPAND_HH */
