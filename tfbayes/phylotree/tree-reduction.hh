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

#ifndef __TFBAYES_PHYLOTREE_TREE_REDUCTION_HH__
#define __TFBAYES_PHYLOTREE_TREE_REDUCTION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>

#include <tfbayes/phylotree/model.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/polynomial.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

namespace tfbayes_detail {
        class leafset_t : public boost::unordered_set<const pt_leaf_t*> {
        public:
                bool empty() const {
                        return size() == 0;
                }
                void join(const leafset_t& leafset) {
                        for (leafset_t::const_iterator it = leafset.begin(); it != leafset.end(); it++) {
                                insert(*it);
                        }
                }
        };
        inline
        size_t hash_value(const leafset_t& set) {
                size_t seed = 0;
                for (leafset_t::const_iterator it = set.begin(); it != set.end(); it++) {
                        seed += (size_t)*it;
                }
                return seed;
        }
        std::ostream& operator<< (std::ostream& o, const leafset_t& leafset) {
                for (leafset_t::const_iterator it = leafset.begin(); it != leafset.end(); it++) {
                        if (it != leafset.begin()) {
                                o << ",";
                        }
                        o << (*it)->name;
                }
                return o;
        }
        class incomplete_leafset_t : public boost::unordered_set<leafset_t> {
        public:
                bool is_complete() const {
                        if (incomplete().empty()) {
                                return true;
                        }
                        return false;
                }
                incomplete_leafset_t& complete() {
                        if (!is_complete()) {
                                insert(incomplete());
                                incomplete() = leafset_t();
                        }
                        return *this;
                }
                leafset_t& incomplete() {
                        return _incomplete;
                }
                const leafset_t& incomplete() const {
                        return _incomplete;
                }
                incomplete_leafset_t& operator*=(const incomplete_leafset_t& leafset) {
                        for (incomplete_leafset_t::const_iterator it = leafset.begin(); it != leafset.end(); it++) {
                                insert(*it);
                        }
                        for (leafset_t::const_iterator it = leafset.incomplete().begin(); it != leafset.incomplete().end(); it++) {
                                incomplete().insert(*it);
                        }
                        return *this;
                }
                bool operator==(const incomplete_leafset_t& leafset) const {
                        return (const boost::unordered_set<leafset_t>&)*this == (const boost::unordered_set<leafset_t>&)leafset &&
                                incomplete() == leafset.incomplete();
                }
        protected:
                leafset_t _incomplete;
        };
        inline
        size_t hash_value(const incomplete_leafset_t& leafset) {
                size_t seed = 0;
                for (incomplete_leafset_t::const_iterator it = leafset.begin(); it != leafset.end(); it++) {
                        seed += hash_value(*it);
                }
                seed += hash_value(leafset.incomplete());
                return seed;
        }
        
        class incomplete_nodeterm_t : public incomplete_leafset_t {
        public:
                incomplete_nodeterm_t()
                        : incomplete_leafset_t(), _coefficient(1.0) { }
                incomplete_nodeterm_t(std::pair<const incomplete_leafset_t, double>& pair)
                        : incomplete_leafset_t(pair.first), _coefficient(pair.second) { }

                incomplete_nodeterm_t& complete() {
                        incomplete_leafset_t::complete();
                        return *this;
                }
                incomplete_nodeterm_t& operator*=(double constant) {
                        coefficient() *= constant;
                        return *this;
                }
                incomplete_nodeterm_t& operator*=(const incomplete_nodeterm_t& term) {
                        incomplete_leafset_t::operator*=(term);
                        coefficient() *= term.coefficient();
                        return *this;
                }
                
                const double& coefficient() const {
                        return _coefficient;
                }
                double& coefficient() {
                        return _coefficient;
                }
                
        private:
                double _coefficient;
        };
        inline
        incomplete_nodeterm_t operator*(double constant, const incomplete_nodeterm_t& term) {
                incomplete_nodeterm_t result(term);
                result *= constant;
                return result;
        }
        inline
        incomplete_nodeterm_t operator*(const incomplete_nodeterm_t& term, double constant) {
                incomplete_nodeterm_t result(term);
                result *= constant;
                return result;
        }
        inline
        incomplete_nodeterm_t operator*(const incomplete_nodeterm_t& term1, const incomplete_nodeterm_t& term2) {
                incomplete_nodeterm_t result(term1);
                result *= term2;
                return result;
        }

        class incomplete_expression_t : public boost::unordered_map<incomplete_leafset_t, double> {
        public:
                std::pair<size_t, size_t> size() {
                        size_t   complete = 0;
                        size_t incomplete = 0;
                        for (incomplete_expression_t::const_iterator it = begin(); it != end(); it++) {
                                if (it->is_complete()) {
                                        complete++;
                                }
                                else {
                                        incomplete++;
                                }
                        }
                        return std::pair<size_t, size_t>(complete, incomplete);
                }
                incomplete_expression_t& operator+=(const incomplete_nodeterm_t& term) {
                        operator[](term) += term.coefficient();
                        if (operator[](term) == 0.0) {
                                erase(term);
                        }
                        return *this;
                }
                incomplete_expression_t& operator-=(const incomplete_nodeterm_t& term) {
                        operator[](term) -= term.coefficient();
                        if (operator[](term) == 0.0) {
                                erase(term);
                        }
                        return *this;
                }
                incomplete_expression_t& operator*=(const incomplete_nodeterm_t& term) {
                        incomplete_expression_t tmp;

                        for (incomplete_expression_t::const_iterator it = this->begin(); it != this->end(); it++) {
                                tmp += (*it)*term;
                        }
                        operator=(tmp);

                        return *this;
                }
                incomplete_expression_t& operator+=(const incomplete_expression_t& expression) {
                        for (incomplete_expression_t::const_iterator it = expression.begin(); it != expression.end(); it++) {
                                operator+=(*it);
                        }
                        return *this;
                }
                incomplete_expression_t& operator-=(const incomplete_expression_t& expression) {
                        for (incomplete_expression_t::const_iterator it = expression.begin(); it != expression.end(); it++) {
                                operator-=(*it);
                        }
                        return *this;
                }
                incomplete_expression_t& operator*=(const incomplete_expression_t& expression) {
                        incomplete_expression_t tmp;

                        for (incomplete_expression_t::const_iterator it = this->begin(); it != this->end(); it++) {
                                for (incomplete_expression_t::const_iterator is = expression.begin(); is != expression.end(); is++) {
                                        tmp += (*it)*(*is);
                                }
                        }
                        operator=(tmp);

                        return *this;
                }

                // Iterator
                ////////////////////////////////////////////////////////////////////////
                class const_iterator : public boost::unordered_map<incomplete_leafset_t, double>::const_iterator
                {
                public:
                        const_iterator(boost::unordered_map<incomplete_leafset_t, double>::const_iterator iterator)
                                : boost::unordered_map<incomplete_leafset_t, double>::const_iterator(iterator)
                                { }

                        const incomplete_nodeterm_t* operator->() const
                                {
                                        return (const incomplete_nodeterm_t*)boost::unordered_map<incomplete_leafset_t, double>::const_iterator::operator->();
                                }
                        const incomplete_nodeterm_t operator*() const
                                {
                                        return *operator->();
                                }
                };
                const_iterator begin() const {
                        return const_iterator(boost::unordered_map<incomplete_leafset_t, double>::begin());
                }
        };
        incomplete_expression_t operator*(const incomplete_expression_t& expression1, const incomplete_expression_t& expression2) {
                incomplete_expression_t result(expression1);
                result *= expression2;
                return result;
        }
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
        inline
        std::ostream& operator<< (std::ostream& o, const incomplete_leafset_t& exponent) {
                for (incomplete_leafset_t::const_iterator it = exponent.begin(); it != exponent.end(); it++) {
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
        inline
        std::ostream& operator<< (std::ostream& o, const incomplete_nodeterm_t& term) {
                if (term.coefficient() != 1.0) {
                        o << term.coefficient()
                          << " ";
                }
                o << (const incomplete_leafset_t&)term;
                return o;
        }
        inline
        std::ostream& operator<< (std::ostream& o, const incomplete_expression_t& expression) {
                for (incomplete_expression_t::const_iterator it = expression.begin(); it != expression.end(); it++) {
                        if (it != expression.begin()) {
                                o << " + ";
                        }
                        o << *it;
                }
                return o;
        }
};

using tfbayes_detail::pt_expand_t;
using tfbayes_detail::incomplete_expression_t;
using tfbayes_detail::operator<<;

incomplete_expression_t pt_simplify(const pt_root_t& node) {
        return tfbayes_detail::pt_simplify_rec(node);
}

#endif /* __TFBAYES_PHYLOTREE_TREE_REDUCTION_HH__ */
