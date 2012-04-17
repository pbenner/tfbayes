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

#ifndef INCOMPLETE_POLYNOMIAL_HH
#define INCOMPLETE_POLYNOMIAL_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <tfbayes/polynomial.hh>

#include <phylotree.hh>
#include <phylotree-expand.hh>
#include <nodeset.hh>

template <typename T, size_t S>
class incomplete_exponent_t : public exponent_t<T, S> {
public:
        incomplete_exponent_t() { }
        incomplete_exponent_t(const exponent_t<T, S>& exponent)
                : exponent_t<T, S>(exponent) { }

        bool is_complete() const {
                if (incomplete().empty()) {
                        return true;
                }
                return false;
        }

        nodeset_t& incomplete() {
                return _incomplete;
        }
        const nodeset_t& incomplete() const {
                return _incomplete;
        }

        incomplete_exponent_t<T, S>& operator*=(const incomplete_exponent_t<T, S>& exponent) {
                exponent_t<T, S>::operator*=(exponent);
                for (nodeset_t::const_iterator it = exponent.incomplete().begin(); it != exponent.incomplete().end(); it++) {
                        incomplete().insert(*it);
                }

                return *this;
        }
        bool operator==(const incomplete_exponent_t<T, S>& exponent) const {
                return static_cast<const exponent_t<T, S>&>(*this) == static_cast<const exponent_t<T, S>&>(exponent) &&
                        incomplete() == exponent.incomplete();
        }

private:
        nodeset_t _incomplete;
};

template <typename T, size_t S> class incomplete_polynomial_t;

template <typename T, size_t S>
class incomplete_term_t : public incomplete_exponent_t<T, S> {
public:
        incomplete_term_t()
                : incomplete_exponent_t<T, S>(), _coefficient(1.0) { }
        incomplete_term_t(std::pair<const incomplete_exponent_t<T, S>, double>& pair)
                : incomplete_exponent_t<T, S>(pair.first), _coefficient(pair.second) { }

        polynomial_t<T, S> complete() const {
                polynomial_t<T, S> poly(1.0);
                polynomial_term_t<T, S> term(coefficient());
                term.exponent() = *this;
                if (!incomplete_exponent_t<T, S>::is_complete()) {
                        poly = pt_expand<T, S>(incomplete_exponent_t<T, S>::incomplete());
                }
                return poly*term;
        }

        incomplete_term_t<T, S>& operator*=(double constant) {
                coefficient() *= constant;
                return *this;
        }
        incomplete_term_t<T, S>& operator*=(const incomplete_term_t<T, S>& term) {
                incomplete_exponent_t<T, S>::operator*=(term);
                coefficient() *= term.coefficient();
                return *this;
        }

        const incomplete_exponent_t<T, S>& exponent() const {
                return *this;
        }
        incomplete_exponent_t<T, S>& exponent() {
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

template <typename T, size_t S>
incomplete_term_t<T, S> operator*(double constant, const incomplete_term_t<T, S>& term) {
        incomplete_term_t<T, S> result(term);
        result *= constant;
        return result;
}
template <typename T, size_t S>
incomplete_term_t<T, S> operator*(const incomplete_term_t<T, S>& term, double constant) {
        incomplete_term_t<T, S> result(term);
        result *= constant;
        return result;
}
template <typename T, size_t S>
incomplete_term_t<T, S> operator*(const incomplete_term_t<T, S>& term1, const incomplete_term_t<T, S>& term2) {
        incomplete_term_t<T, S> result(term1);
        result *= term2;
        return result;
}

template <typename T, size_t S>
class incomplete_polynomial_t : public boost::unordered_map<incomplete_exponent_t<T, S>, double> {
public:
        incomplete_polynomial_t() { }
        incomplete_polynomial_t(const polynomial_t<T, S>& poly) {
                for (typename polynomial_t<T, S>::const_iterator it = poly.begin(); it != poly.end(); it++) {
                        incomplete_term_t<T, S> term;
                        term.coefficient() = (*it).coefficient();
                        term.exponent()    = (*it).exponent();
                        operator+=(term);
                }
        }
        std::pair<size_t, size_t> size() {
                size_t   complete = 0;
                size_t incomplete = 0;
                for (typename incomplete_polynomial_t<T, S>::const_iterator it = this->begin(); it != this->end(); it++) {
                        if (it->is_complete()) {
                                complete++;
                        }
                        else {
                                incomplete++;
                        }
                }
                return std::pair<size_t, size_t>(complete, incomplete);
        }
        polynomial_t<T, S> complete() {
                polynomial_t<T, S> poly;
                for (typename incomplete_polynomial_t<T, S>::const_iterator it = this->begin(); it != this->end(); it++) {
                        poly += (*it).complete();
                }
                return poly;
        }
        incomplete_polynomial_t& operator+=(const incomplete_term_t<T, S>& term) {
                operator[](term) += term.coefficient();
                if (operator[](term) == 0.0) {
                        erase(term);
                }
                return *this;
        }
        incomplete_polynomial_t& operator-=(const incomplete_term_t<T, S>& term) {
                operator[](term) -= term.coefficient();
                if (operator[](term) == 0.0) {
                        erase(term);
                }
                return *this;
        }
        incomplete_polynomial_t& operator*=(const incomplete_term_t<T, S>& term) {
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
        class const_iterator : public boost::unordered_map<incomplete_exponent_t<T, S>, double>::const_iterator
        {
        public:
                const_iterator(typename boost::unordered_map<incomplete_exponent_t<T, S>, double>::const_iterator iterator)
                        : boost::unordered_map<incomplete_exponent_t<T, S>, double>::const_iterator(iterator)
                        { }

                const incomplete_term_t<T, S>* operator->() const
                {
                        return (const incomplete_term_t<T, S>*)boost::unordered_map<incomplete_exponent_t<T, S>, double>::const_iterator::operator->();
                }
                const incomplete_term_t<T, S> operator*() const
                {
                        return *operator->();
                }
        };
        const_iterator begin() const {
                return const_iterator(boost::unordered_map<incomplete_exponent_t<T, S>, double>::begin());
        }
};

template <typename T, size_t S>
incomplete_polynomial_t<T, S> operator*(
        const incomplete_polynomial_t<T, S>& poly1,
        const incomplete_polynomial_t<T, S>& poly2)
{
        incomplete_polynomial_t<T, S> result(poly1);
        result *= poly2;
        return result;
}

////////////////////////////////////////////////////////////////////////////////

#include <ostream>

template <typename T, size_t S>
std::ostream& operator<< (std::ostream& o, const incomplete_exponent_t<T, S>& exponent) {
        o << (const exponent_t<T, S>)exponent;
        if (!exponent.incomplete().empty()) {
                o << "I("
                  << exponent.incomplete()
                  << ")";
        }
        return o;
}

template <typename T, size_t S>
std::ostream& operator<< (std::ostream& o, const incomplete_term_t<T, S>& term) {
        if (term.coefficient() != 1.0) {
                o << term.coefficient()
                  << " ";
        }
        o << (const incomplete_exponent_t<T, S>&)term;

        return o;
}

template <typename T, size_t S>
std::ostream& operator<< (std::ostream& o, const incomplete_polynomial_t<T, S>& polynomial) {
        for (typename incomplete_polynomial_t<T, S>::const_iterator it = polynomial.begin(); it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + ";
                }
                o << *it;
        }

        return o;
}

#include <boost/functional/hash.hpp> 

template <typename T, size_t S>
size_t hash_value(const incomplete_exponent_t<T, S>& exponent) {
        size_t seed = 0;
        boost::hash_combine(seed, static_cast<const exponent_t<T, S>&>(exponent));
        boost::hash_combine(seed, hash_value(exponent.incomplete()));
        return seed;
}

#endif /* INCOMPLETE_POLYNOMIAL_HH */
