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

#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <utility>

#include <boost/array.hpp>
#include <boost/unordered_map.hpp>

template <typename T, size_t S>
class exponent_t : public boost::array<T, S> {
public:
        exponent_t() {
                for (size_t i = 0; i < S; i++) {
                        operator[](i) = 0;
                }
        }

        exponent_t<T, S>& operator*=(const exponent_t<T, S>& exponent) {
                for (size_t i = 0; i < S; i++) {
                        operator[](i) += exponent[i];
                }
                return *this;
        }
        exponent_t<T, S>& operator/=(const exponent_t<T, S>& exponent) {
                for (size_t i = 0; i < S; i++) {
                        operator[](i) -= exponent[i];
                }
                return *this;
        }
        double eval(const boost::array<double, S>& val) const {
                double result = 1.0;
                for (size_t i = 0; i < S; i++) {
                        result *= pow(val[i], operator[](i));
                }
                return result;
        }
        using boost::array<T, S>::operator[];
};

template <typename T, size_t S, typename C = double>
class polynomial_term_t : public std::pair<exponent_t<T, S>, C> {
public:
        polynomial_term_t()
                : std::pair<exponent_t<T, S>, C>(exponent_t<T, S>(), 0.0)
                { }
        polynomial_term_t(const C& constant)
                : std::pair<exponent_t<T, S>, C>(exponent_t<T, S>(), constant)
                { }
        polynomial_term_t(const std::pair<exponent_t<T, S>, C>& pair)
                : std::pair<exponent_t<T, S>, C>(pair)
                { }
        exponent_t<T, S>& exponent() {
                return std::pair<exponent_t<T, S>, C>::first;
        }
        const exponent_t<T, S>& exponent() const {
                return std::pair<exponent_t<T, S>, C>::first;
        }
        C& coefficient() {
                return std::pair<exponent_t<T, S>, C>::second;
        }
        const C& coefficient() const {
                return std::pair<exponent_t<T, S>, C>::second;
        }

        polynomial_term_t<T, S, C> operator-() const {
                return std::pair<exponent_t<T, S>, C>(exponent(), -coefficient());
        }
        polynomial_term_t<T, S, C>& operator+=(const polynomial_term_t<T, S, C>& term) {
                coefficient() += term.coefficient();
                return *this;
        }
        polynomial_term_t<T, S, C>& operator-=(const polynomial_term_t<T, S, C>& term) {
                coefficient() -= term.coefficient();
                return *this;
        }
        polynomial_term_t<T, S, C>& operator*=(const C& constant) {
                coefficient() *= constant;
                return *this;
        }
        polynomial_term_t<T, S, C>& operator/=(const C& constant) {
                coefficient() /= constant;
                return *this;
        }
        polynomial_term_t<T, S, C>& operator*=(const polynomial_term_t<T, S, C>& term) {
                coefficient() *= term.coefficient();
                exponent()    *= term.exponent();
                return *this;
        }
        polynomial_term_t<T, S, C>& operator/=(const polynomial_term_t<T, S, C>& term) {
                coefficient() /= term.coefficient();
                exponent()    /= term.exponent();
                return *this;
        }
        C eval(const boost::array<double, S>& val) const {
                C c(coefficient());
                c *= exponent().eval(val);
                return c;
        }
};

template <typename T, size_t S, typename C>
polynomial_term_t<T, S, C> operator+(const polynomial_term_t<T, S, C>& term1, const polynomial_term_t<T, S, C>& term2) {
        polynomial_term_t<T, S, C> tmp(term1);
        tmp += term2;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_term_t<T, S, C> operator-(const polynomial_term_t<T, S, C>& term1, const polynomial_term_t<T, S, C>& term2) {
        polynomial_term_t<T, S, C> tmp(term1);
        tmp -= term2;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_term_t<T, S, C> operator*(const polynomial_term_t<T, S, C>& term, const C& constant) {
        polynomial_term_t<T, S, C> tmp(term);
        tmp *= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_term_t<T, S, C> operator/(const polynomial_term_t<T, S, C>& term, const C& constant) {
        polynomial_term_t<T, S, C> tmp(term);
        tmp /= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_term_t<T, S, C> operator*(const C& constant, const polynomial_term_t<T, S, C>& term) {
        polynomial_term_t<T, S, C> tmp(term);
        tmp *= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_term_t<T, S, C> operator/(const C& constant, const polynomial_term_t<T, S, C>& term) {
        polynomial_term_t<T, S, C> tmp(constant);
        tmp /= term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_term_t<T, S, C> operator*(const polynomial_term_t<T, S, C>& term1, const polynomial_term_t<T, S, C>& term2) {
        polynomial_term_t<T, S, C> tmp(term1);
        tmp *= term2;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_term_t<T, S, C> operator/(const polynomial_term_t<T, S, C>& term1, const polynomial_term_t<T, S, C>& term2) {
        polynomial_term_t<T, S, C> tmp(term1);
        tmp /= term2;
        return tmp;
}

template <typename T, size_t S, typename C = double>
class polynomial_t : public boost::unordered_map<exponent_t<T, S>, C> {
public:
        polynomial_t()
                : boost::unordered_map<exponent_t<T, S>, C>() {}
        polynomial_t(const C& constant)
                : boost::unordered_map<exponent_t<T, S>, C>() {
                operator+=(polynomial_term_t<T, S, C>(constant));
        }
        polynomial_t(const polynomial_term_t<T, S, C>& term)
                : boost::unordered_map<exponent_t<T, S>, C>() {
                operator+=(term);
        }

        polynomial_t<T, S, C>& operator+=(const C& constant) {
                operator+=(polynomial_term_t<T, S, C>(constant));
                return *this;
        }
        polynomial_t<T, S, C>& operator-=(const C& constant) {
                operator-=(polynomial_term_t<T, S, C>(constant));
                return *this;
        }
        polynomial_t<T, S, C>& operator+=(const polynomial_term_t<T, S, C>& term) {
                if (term.coefficient()) {
                        operator[](term.exponent()) += term.coefficient();
                        if (operator[](term.exponent()) == 0.0) {
                                erase(term.exponent());
                        }
                }
                return *this;
        }
        polynomial_t<T, S, C>& operator-=(const polynomial_term_t<T, S, C>& term) {
                if (term.coefficient() != 0.0) {
                        operator[](term.exponent()) -= term.coefficient();
                        if (operator[](term.exponent()) == 0.0) {
                                erase(term.exponent());
                        }
                }
                return *this;
        }
        polynomial_t<T, S, C>& operator+=(const polynomial_t<T, S, C>& poly) {
                for (const_iterator it = poly.begin(); it != poly.end(); it++) {
                        operator+=(*it);
                }
                return *this;
        }
        polynomial_t<T, S, C>& operator-=(const polynomial_t<T, S, C>& poly) {
                for (const_iterator it = poly.begin(); it != poly.end(); it++) {
                        operator-=(*it);
                }
                return *this;
        }
        polynomial_t<T, S, C>& operator*=(const C& constant) {
                polynomial_t<T, S, C> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) * constant;
                }
                operator=(tmp);

                return *this;
        }
        polynomial_t<T, S, C>& operator/=(const C& constant) {
                polynomial_t<T, S, C> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) / constant;
                }
                operator=(tmp);

                return *this;
        }
        polynomial_t<T, S, C>& operator*=(const polynomial_term_t<T, S, C>& term) {
                polynomial_t<T, S, C> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) * term;
                }
                operator=(tmp);

                return *this;
        }
        polynomial_t<T, S, C>& operator/=(const polynomial_term_t<T, S, C>& term) {
                polynomial_t<T, S, C> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) / term;
                }
                operator=(tmp);

                return *this;
        }
        polynomial_t<T, S, C>& operator*=(const polynomial_t<T, S, C>& poly) {
                polynomial_t<T, S, C> tmp;

                for (typename polynomial_t<T, S, C>::const_iterator it = this->begin(); it != this->end(); it++) {
                        for (typename polynomial_t<T, S, C>::const_iterator is = poly.begin(); is != poly.end(); is++) {
                                tmp += (*it)*(*is);
                        }
                }
                operator=(tmp);

                return *this;
        }
        C eval(const boost::array<double, S>& val) const {
                C result(0.0);
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        result += it->eval(val);
                }
                return result;
        }
        using boost::unordered_map<exponent_t<T, S>, C>::operator[];
        using boost::unordered_map<exponent_t<T, S>, C>::erase;

        // Iterator
        ////////////////////////////////////////////////////////////////////////
        class const_iterator : public boost::unordered_map<exponent_t<T, S>, C>::const_iterator
        {
        public:
                const_iterator(typename boost::unordered_map<exponent_t<T, S>, C>::const_iterator iterator)
                        : boost::unordered_map<exponent_t<T, S>, C>::const_iterator(iterator)
                        {}

                const polynomial_term_t<T, S, C>* operator->() const
                {
                        return (const polynomial_term_t<T, S, C>*)boost::unordered_map<exponent_t<T, S>, C>::const_iterator::operator->();
                }
                const polynomial_term_t<T, S, C> operator*() const
                {
                        return *operator->();
                }
        };
        const_iterator begin() const {
                return const_iterator(boost::unordered_map<exponent_t<T, S>, C>::begin());
        }
};

template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator+(const polynomial_t<T, S, C>& poly, const C& constant) {
        polynomial_t<T, S, C> tmp(poly);
        tmp += constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator-(const polynomial_t<T, S, C>& poly, const C& constant) {
        polynomial_t<T, S, C> tmp(poly);
        tmp -= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator+(const C& constant, const polynomial_t<T, S, C>& poly) {
        polynomial_t<T, S, C> tmp(poly);
        tmp += constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator-(const C& constant, const polynomial_t<T, S, C>& poly) {
        polynomial_t<T, S, C> tmp(poly);
        tmp -= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator*(const polynomial_t<T, S, C>& poly, const C& constant) {
        polynomial_t<T, S, C> tmp(poly);
        tmp *= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator/(const polynomial_t<T, S, C>& poly, const C& constant) {
        polynomial_t<T, S, C> tmp(poly);
        tmp /= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator*(const C& constant, const polynomial_t<T, S, C>& poly) {
        polynomial_t<T, S, C> tmp(poly);
        tmp *= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator/(const C& constant, const polynomial_t<T, S, C>& poly) {
        polynomial_t<T, S, C> tmp(poly);
        tmp /= constant;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator+(const polynomial_t<T, S, C>& poly, const polynomial_term_t<T, S, C>& term) {
        polynomial_t<T, S, C> tmp(poly);
        tmp += term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator-(const polynomial_t<T, S, C>& poly, const polynomial_term_t<T, S, C>& term) {
        polynomial_t<T, S, C> tmp(poly);
        tmp -= term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator+(const polynomial_term_t<T, S, C>& term, const polynomial_t<T, S, C>& poly) {
        polynomial_t<T, S, C> tmp(poly);
        tmp += term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator-(const polynomial_term_t<T, S, C>& term, const polynomial_t<T, S, C>& poly) {
        polynomial_t<T, S, C> tmp(poly);
        tmp -= term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator*(const polynomial_t<T, S, C>& poly, const polynomial_term_t<T, S, C>& term) {
        polynomial_t<T, S, C> tmp(poly);
        tmp *= term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator/(const polynomial_t<T, S, C>& poly, const polynomial_term_t<T, S, C>& term) {
        polynomial_t<T, S, C> tmp(poly);
        tmp /= term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator*(const polynomial_term_t<T, S, C>& term, const polynomial_t<T, S, C>& poly) {
        polynomial_t<T, S, C> tmp(poly);
        tmp *= term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator/(const polynomial_term_t<T, S, C>& term, const polynomial_t<T, S, C>& poly) {
        polynomial_t<T, S, C> tmp(poly);
        tmp /= term;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator+(const polynomial_t<T, S, C>& poly1, const polynomial_t<T, S, C>& poly2) {
        polynomial_t<T, S, C> tmp(poly1);
        tmp += poly2;
        return tmp;
}
template <typename T, size_t S, typename C>
polynomial_t<T, S, C> operator*(const polynomial_t<T, S, C>& poly1, const polynomial_t<T, S, C>& poly2) {
        if (poly1.size() == 0) {
                return polynomial_t<T, S, C>();
        }
        if (poly2.size() == 0) {
                return polynomial_t<T, S, C>();
        }
        polynomial_t<T, S, C> tmp(poly1);
        tmp *= poly2;
        return tmp;
}

#endif /* _POLYNOMIAL_H_ */
