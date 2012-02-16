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
        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> px;
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
polynomial_t<CODE_TYPE, ALPHABET_SIZE> pt_expand(const node_set_t& node_set) {
        std::vector<bool> applicable(ALPHABET_SIZE, false);
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> result;
        polynomial_t<CODE_TYPE, ALPHABET_SIZE> remainder(1.0);

        for (node_set_t::const_iterator it = node_set.begin(); it != node_set.end(); it++) {
                applicable[(*it)->x] = true;
        }
        for (CODE_TYPE y = 0; y < ALPHABET_SIZE; y++) {
                if (applicable[y]) {
                        polynomial_term_t<CODE_TYPE, ALPHABET_SIZE> px = nucleotide_probability<CODE_TYPE, ALPHABET_SIZE>(y);
                        polynomial_t<CODE_TYPE, ALPHABET_SIZE> tmp(px);
                        remainder -= px;
                        for (node_set_t::const_iterator it = node_set.begin(); it != node_set.end(); it++) {
                                tmp *= mutation_model<CODE_TYPE, ALPHABET_SIZE>(*it, y);
                        }
                        result += tmp;
                }
        }
        if (remainder.size() < ALPHABET_SIZE) {
                polynomial_t<CODE_TYPE, ALPHABET_SIZE> tmp(remainder);
                for (node_set_t::const_iterator it = node_set.begin(); it != node_set.end(); it++) {
                        tmp *= mutation_model<CODE_TYPE, ALPHABET_SIZE>(*it, ALPHABET_SIZE);
                }
                result += tmp;
        }
        return result;
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
