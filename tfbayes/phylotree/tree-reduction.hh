/* Copyright (C) 2012-2013 Philipp Benner
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

#ifndef __TFBAYES_PHYLOTREE_PHYLOTREE_EXPAND_HH__
#define __TFBAYES_PHYLOTREE_PHYLOTREE_EXPAND_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/incomplete-expression.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/polynomial.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

incomplete_expression_t pt_simplify_leaf(const pt_leaf_t& node) {
        incomplete_expression_t expression;
        incomplete_nodeterm_t term;
        term.incomplete().insert(&node);
        expression += term;

        return expression;
}

incomplete_expression_t pt_simplify_root(
        const pt_node_t& node,
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
        const pt_node_t& node,
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
                        expression2 += (1.0-node.mutation_probability())*term;
                        expression2 += node.mutation_probability() *term.complete();
                }
        }
        return expression2;
}

incomplete_expression_t pt_simplify_rec(const pt_node_t& node) {
        if (node.leaf()) {
                return pt_simplify_leaf(static_cast<const pt_leaf_t&>(node));
        }
        else {
                const incomplete_expression_t expression_left = pt_simplify_rec(node.left ());
                const incomplete_expression_t expression_right = pt_simplify_rec(node.right());

                if (node.root()) {
                        return pt_simplify_root(node, expression_left, expression_right);
                }
                else {
                        return pt_simplify_node(node, expression_left, expression_right);
                }
        }
}

incomplete_expression_t pt_simplify(const pt_root_t& node) {
        return pt_simplify_rec(node);
}

template <size_t AS, typename AC, typename PC>
polynomial_term_t<AS, PC> nucleotide_probability(AC x) {
        polynomial_term_t<AS, PC> px(1.0);
        px.exponent()[x] = 1;
        return px;
}

template <size_t AS, typename AC, typename PC>
polynomial_t<AS, PC> mutation_model(const pt_node_t& node, AC x, AC y) {
        polynomial_t<AS, PC> poly;
        polynomial_term_t<AS, PC> py = nucleotide_probability<AS, AC, PC>(y);
        poly += node.mutation_probability()*py;
        if (x == y) {
                poly += (1.0-node.mutation_probability());
        }
        return poly;
}

template <size_t AS, typename AC = alphabet_code_t, typename PC = double>
class pt_expand_t : public polynomial_t<AS, PC> {
public:
        using polynomial_t<AS, PC>::operator=;

        pt_expand_t(
                const incomplete_expression_t& expression,
                const std::vector<AC>& observations) {

                polynomial_t<AS, PC> result;
                for (incomplete_expression_t::const_iterator it = expression.begin();
                     it != expression.end(); it++) {
                        result += pt_expand(*it, observations);
                }
                
                operator=(result);
        }
        pt_expand_t(
                const incomplete_nodeterm_t& nodeterm,
                const std::vector<AC>& observations) {

                operator=(pt_expand(nodeterm, observations));
        }
private:
        /*
         * phi(y  ; x) = P(nM) : if x == y; 0 : otherwise
         * phi(nil; x) = P( M) P(x)
         *
         * phi(y  ; x_1, x_2, ..., x_n) = 1[x_n == y] P(nM) { phi(y  ; x_1, x_2, ..., x_{n-1})   +
         *                                                    phi(nil; x_1, x_2, ..., x_{n-1}) } +
         *                                P(M) P(x_n) phi(y  ; x_1, x_2, ..., x_{n-1})
         * phi(nil; x_1, x_2, ..., x_n) = P(M) P(x_n) phi(nil; x_1, x_2, ..., x_{n-1})
         */
        static
        polynomial_t<AS, PC> pt_expand_rec(
                leafset_t::const_iterator it,
                leafset_t::const_iterator end,
                const std::vector<AC>& observations,
                AC condition) {

                polynomial_t<AS, PC> result(0.0);

                if(it == end) {
                        if (condition == AS) {
                                result += 1.0;
                        }
                }
                else {
                        const AC x = observations[(*it)->id];
                        const double pm   = (*it)->mutation_probability();
                        const polynomial_term_t<AS, PC> px =
                                nucleotide_probability<AS, AC, PC>(x);
                        
                        polynomial_t<AS, PC> tmp = pt_expand_rec(++it, end, observations, condition);
                        
                        result += pm*px*tmp;
                        
                        if (condition == x) {
                                result += (1.0-pm)*tmp;
                                result += (1.0-pm)*pt_expand_rec(it, end, observations, AS);
                        }
                }
                return result;
        }
        
        /*
         * phi(x_1, x_2, ..., x_n) = phi(nil; x_1, x_2, x_{n-1})
         *                         + sum_y p(y) phi(y; x_1, x_2, ..., x_{n-1})
         */
        static
        polynomial_t<AS, PC> pt_expand_rec(
                leafset_t::const_iterator it,
                leafset_t::const_iterator end,
                const std::vector<AC>& observations) {

                polynomial_t<AS, PC> result(0.0);
                        
                for (size_t x = 0; x < AS; x++) {
                        const polynomial_term_t<AS, PC> px =
                                nucleotide_probability<AS, AC, PC>(x);
                        result += px*pt_expand_rec(it, end, observations, x);
                }
                result += pt_expand_rec(it, end, observations, AS);
                
                return result;
        }
        
        static
        polynomial_t<AS, PC> pt_expand_rec(
                polynomial_term_t<AS, PC> term,
                leafset_t::const_iterator it,
                leafset_t::const_iterator end,
                const std::vector<AC>& observations,
                AC condition) {

                polynomial_t<AS, PC> result(0.0);
                if(it == end) {
                        /* leaf */
                        if (condition != AS) {
                                term *= nucleotide_probability<AS, AC, PC>(condition);
                        }
                        result += term;
                }
                else {
                        const double pm   = (*it)->mutation_probability();
                        const AC x = observations[(*it)->id]; it++;
                        
                        /* no mutation */
                        if (condition == AS || condition == x) {
                                result += pt_expand_rec((1.0-pm)*term, it, end, x);
                        }
                        /* mutation */
                        term   *= nucleotide_probability<AS, AC, PC>(x);
                        result += pt_expand_rec(pm*term, it, end, condition);
                }
                return result;
        }
        
        static
        polynomial_t<AS, PC> pt_expand(
                const leafset_t& leafset,
                const std::vector<AC>& observations) {

                /* Algorithm 1 */
                // return pt_expand_rec(1.0, leafset.begin(), leafset.end(), AS);
                /* Algorithm 2 */
                return pt_expand_rec(leafset.begin(), leafset.end(), observations);
        }
        static
        polynomial_t<AS, PC> pt_expand(
                const incomplete_nodeterm_t& nodeterm,
                const std::vector<AC>& observations) {

                polynomial_t<AS, PC> result(1.0);
                for (incomplete_nodeterm_t::const_iterator it = nodeterm.begin();
                     it != nodeterm.end(); it++) {
                        result *= pt_expand(*it, observations);
                }
                return nodeterm.coefficient()*result;
        }
};

#endif /* __TFBAYES_PHYLOTREE_PHYLOTREE_EXPAND_HH__ */
