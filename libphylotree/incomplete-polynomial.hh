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

#ifndef INCOMPLETE_POLYNOMIAL_H
#define INCOMPLETE_POLYNOMIAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <phylotree.hh>
#include <utility.hh>

class node_set_t : public boost::unordered_set<pt_node_t<code_t, alphabet_size>*> {
public:
        bool empty() const {
                return size() == 0;
        }
};

class incomplete_exponent_t : public boost::unordered_set<node_set_t> {
public:
        incomplete_exponent_t& complete() {
                insert(incomplete());
                incomplete() = node_set_t();
                return *this;
        }

        node_set_t& incomplete() {
                return _incomplete;
        }
        const node_set_t& incomplete() const {
                return _incomplete;
        }

        incomplete_exponent_t& operator*=(const incomplete_exponent_t& exponent) {
                for (incomplete_exponent_t::const_iterator it = exponent.begin(); it != exponent.end(); it++) {
                        insert(*it);
                }
                for (node_set_t::const_iterator it = exponent.incomplete().begin(); it != exponent.incomplete().end(); it++) {
                        incomplete().insert(*it);
                }

                return *this;
        }
        bool operator==(const incomplete_exponent_t& exponent) const {
                return (const boost::unordered_set<node_set_t>&)*this == (const boost::unordered_set<node_set_t>&)exponent &&
                       incomplete() == exponent.incomplete();
        }

private:
        node_set_t _incomplete;
};

class incomplete_term_t : public incomplete_exponent_t {
public:
        incomplete_term_t()
                : incomplete_exponent_t(), _coefficient(1.0) { }
        incomplete_term_t(std::pair<const incomplete_exponent_t, double>& pair)
                : incomplete_exponent_t(pair.first), _coefficient(pair.second) { }

        incomplete_term_t& complete() {
                incomplete_exponent_t::complete();
                return *this;
        }

        incomplete_term_t& operator*=(double constant) {
                coefficient() *= constant;
                return *this;
        }
        incomplete_term_t& operator*=(const incomplete_term_t& term) {
                incomplete_exponent_t::operator*=(term);
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

incomplete_term_t operator*(double constant, const incomplete_term_t& term) {
        incomplete_term_t result(term);
        result *= constant;
        return result;
}
incomplete_term_t operator*(const incomplete_term_t& term, double constant) {
        incomplete_term_t result(term);
        result *= constant;
        return result;
}
incomplete_term_t operator*(const incomplete_term_t& term1, const incomplete_term_t& term2) {
        incomplete_term_t result(term1);
        result *= term2;
        return result;
}

class incomplete_polynomial_t : public boost::unordered_map<incomplete_exponent_t, double> {
public:
        incomplete_polynomial_t& operator+=(const incomplete_term_t& term) {
                operator[](term) += term.coefficient();
                if (operator[](term) == 0.0) {
                        erase(term);
                }
                return *this;
        }
        incomplete_polynomial_t& operator-=(const incomplete_term_t& term) {
                operator[](term) -= term.coefficient();
                if (operator[](term) == 0.0) {
                        erase(term);
                }
                return *this;
        }
        incomplete_polynomial_t& operator*=(const incomplete_term_t& term) {
                incomplete_polynomial_t tmp;

                for (incomplete_polynomial_t::const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it)*term;
                }
                operator=(tmp);
        
                return *this;
        }
        incomplete_polynomial_t& operator+=(const incomplete_polynomial_t& poly) {
                for (incomplete_polynomial_t::const_iterator it = poly.begin(); it != poly.end(); it++) {
                        operator+=(*it);
                }
                return *this;
        }
        incomplete_polynomial_t& operator-=(const incomplete_polynomial_t& poly) {
                for (incomplete_polynomial_t::const_iterator it = poly.begin(); it != poly.end(); it++) {
                        operator-=(*it);
                }
                return *this;
        }
        incomplete_polynomial_t& operator*=(const incomplete_polynomial_t& poly) {
                incomplete_polynomial_t tmp;

                for (incomplete_polynomial_t::const_iterator it = this->begin(); it != this->end(); it++) {
                        for (incomplete_polynomial_t::const_iterator is = poly.begin(); is != poly.end(); is++) {
                                tmp += (*it)*(*is);
                        }
                }
                operator=(tmp);
        
                return *this;
        }

        // Iterator
        ////////////////////////////////////////////////////////////////////////
        class const_iterator : public boost::unordered_map<incomplete_exponent_t, double>::const_iterator
        {
        public:
                const_iterator(boost::unordered_map<incomplete_exponent_t, double>::const_iterator iterator)
                        : boost::unordered_map<incomplete_exponent_t, double>::const_iterator(iterator)
                        { }

                const incomplete_term_t* operator->() const
                {
                        return (const incomplete_term_t*)boost::unordered_map<incomplete_exponent_t, double>::const_iterator::operator->();
                }
                const incomplete_term_t operator*() const
                {
                        return *operator->();
                }
        };
        const_iterator begin() const {
                return const_iterator(boost::unordered_map<incomplete_exponent_t, double>::begin());
        }
};

incomplete_polynomial_t operator*(const incomplete_polynomial_t& poly1, const incomplete_polynomial_t& poly2) {
        incomplete_polynomial_t result(poly1);
        result *= poly2;
        return result;
}

#endif /* INCOMPLETE_POLYNOMIAL_H */
