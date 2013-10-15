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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/incomplete-expression.hh>
#include <tfbayes/utility/polynomial.hh>

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_term_t<ALPHABET_SIZE, CODE_TYPE> nucleotide_probability(CODE_TYPE x) {
        polynomial_term_t<ALPHABET_SIZE, CODE_TYPE> px(1.0);
        px.exponent()[x] = 1;
        return px;
}

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE> mutation_model(const pt_node_t* node, CODE_TYPE x, CODE_TYPE y) {
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> poly;
        polynomial_term_t<ALPHABET_SIZE, CODE_TYPE> py = nucleotide_probability<ALPHABET_SIZE, CODE_TYPE>(y);
        poly += node->mutation_probability()*py;
        if (x == y) {
                poly += (1.0-node->mutation_probability());
        }
        return poly;
}

/*
 * phi(y  ; x) = P(nM) : if x == y; 0 : otherwise
 * phi(nil; x) = P( M) P(x)
 *
 * phi(y  ; x_1, x_2, ..., x_n) = 1[x_n == y] P(nM) { phi(y  ; x_1, x_2, ..., x_{n-1})   +
 *                                                    phi(nil; x_1, x_2, ..., x_{n-1}) } +
 *                                P(M) P(x_n) phi(y  ; x_1, x_2, ..., x_{n-1})
 * phi(nil; x_1, x_2, ..., x_n) = P(M) P(x_n) phi(nil; x_1, x_2, ..., x_{n-1})
 */
template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE> pt_expand_rec(
        leafset_t::const_iterator it,
        leafset_t::const_iterator end,
        const std::vector<CODE_TYPE>& observations,
        CODE_TYPE condition)
{
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> result(0.0);

        if(it == end) {
                if (condition == ALPHABET_SIZE) {
                        result += 1.0;
                }
        }
        else {
                const CODE_TYPE x = observations[(*it)->id];
                const double pm   = (*it)->mutation_probability();
                const polynomial_term_t<ALPHABET_SIZE, CODE_TYPE> px =
                        nucleotide_probability<ALPHABET_SIZE, CODE_TYPE>(x);

                polynomial_t<ALPHABET_SIZE, CODE_TYPE> tmp = pt_expand_rec<ALPHABET_SIZE, CODE_TYPE>(++it, end, observations, condition);

                result += pm*px*tmp;

                if (condition == x) {
                        result += (1.0-pm)*tmp;
                        result += (1.0-pm)*pt_expand_rec<ALPHABET_SIZE, CODE_TYPE>(it, end, observations, ALPHABET_SIZE);
                }
        }
        return result;
}

/*
 * phi(x_1, x_2, ..., x_n) = phi(nil; x_1, x_2, x_{n-1})
 *                         + sum_y p(y) phi(y; x_1, x_2, ..., x_{n-1})
 */
template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE> pt_expand_rec(
        leafset_t::const_iterator it,
        leafset_t::const_iterator end,
        const std::vector<CODE_TYPE>& observations)
{
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> result(0.0);

        for (size_t x = 0; x < ALPHABET_SIZE; x++) {
                const polynomial_term_t<ALPHABET_SIZE, CODE_TYPE> px =
                        nucleotide_probability<ALPHABET_SIZE, CODE_TYPE>(x);
                result += px*pt_expand_rec<ALPHABET_SIZE, CODE_TYPE>(it, end, observations, x);
        }
        result += pt_expand_rec<ALPHABET_SIZE, CODE_TYPE>(it, end, observations, ALPHABET_SIZE);

        return result;
}

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE> pt_expand_rec(
        polynomial_term_t<ALPHABET_SIZE, CODE_TYPE> term,
        leafset_t::const_iterator it,
        leafset_t::const_iterator end,
        const std::vector<CODE_TYPE>& observations,
        CODE_TYPE condition)
{
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> result(0.0);
        if(it == end) {
                /* leaf */
                if (condition != ALPHABET_SIZE) {
                        term *= nucleotide_probability<ALPHABET_SIZE, CODE_TYPE>(condition);
                }
                result += term;
        }
        else {
                const double pm   = (*it)->mutation_probability();
                const CODE_TYPE x = observations[(*it)->id]; it++;

                /* no mutation */
                if (condition == ALPHABET_SIZE || condition == x) {
                        result += pt_expand_rec((1.0-pm)*term, it, end, x);
                }
                /* mutation */
                term   *= nucleotide_probability<ALPHABET_SIZE, CODE_TYPE>(x);
                result += pt_expand_rec(pm*term, it, end, condition);
        }
        return result;
}

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE> pt_expand(
        const leafset_t& leafset,
        const std::vector<CODE_TYPE>& observations)
{
        /* Algorithm 1 */
//        return pt_expand_rec<ALPHABET_SIZE, CODE_TYPE>(1.0, leafset.begin(), leafset.end(), ALPHABET_SIZE);
        /* Algorithm 2 */
        return pt_expand_rec<ALPHABET_SIZE, CODE_TYPE>(leafset.begin(), leafset.end(), observations);
}

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE> pt_expand(
        const incomplete_nodeterm_t& nodeterm,
        const std::vector<CODE_TYPE>& observations)
{
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> result(1.0);
        for (incomplete_nodeterm_t::const_iterator it = nodeterm.begin(); it != nodeterm.end(); it++) {
                result *= pt_expand<ALPHABET_SIZE, CODE_TYPE>(*it, observations);
        }
        return nodeterm.coefficient()*result;
}

template <size_t ALPHABET_SIZE, typename CODE_TYPE>
polynomial_t<ALPHABET_SIZE, CODE_TYPE> pt_expand(
        const incomplete_expression_t& expression,
        const std::vector<CODE_TYPE>& observations)
{
        polynomial_t<ALPHABET_SIZE, CODE_TYPE> result;
        for (incomplete_expression_t::const_iterator it = expression.begin(); it != expression.end(); it++) {
                result += pt_expand<ALPHABET_SIZE, CODE_TYPE>(*it, observations);
        }

        return result;
}

#endif /* PHYLOTREE_EXPAND_HH */
