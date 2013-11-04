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

#ifndef __TFBAYES_PHYLOTREE_INCOMPLETE_EXPRESSION_HH__
#define __TFBAYES_PHYLOTREE_INCOMPLETE_EXPRESSION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/leafset.hh>

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

private:
        leafset_t _incomplete;
};

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

incomplete_nodeterm_t operator*(double constant, const incomplete_nodeterm_t& term);
incomplete_nodeterm_t operator*(const incomplete_nodeterm_t& term, double constant);
incomplete_nodeterm_t operator*(const incomplete_nodeterm_t& term1, const incomplete_nodeterm_t& term2);

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

incomplete_expression_t operator*(const incomplete_expression_t& expression1, const incomplete_expression_t& expression2);

////////////////////////////////////////////////////////////////////////////////

#include <ostream>

std::ostream& operator<< (std::ostream& o, const incomplete_nodeterm_t& term);
std::ostream& operator<< (std::ostream& o, const incomplete_expression_t& expression);

size_t hash_value(const incomplete_leafset_t& leafset);

#endif /* __TFBAYES_PHYLOTREE_INCOMPLETE_EXPRESSION_HH__ */
