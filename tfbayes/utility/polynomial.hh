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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/utility/logarithmetic.h>
#include <utility>

#include <boost/array.hpp>
#include <boost/unordered_map.hpp>

template <size_t S, typename T>
class exponent_t : public boost::array<T, S> {
public:
        exponent_t() {
                for (size_t i = 0; i < S; i++) {
                        operator[](i) = 0;
                }
        }
        template <typename InputIterator>
        exponent_t(InputIterator first, InputIterator last)
                : boost::array<T, S>() {
                for (size_t i = 0; i < S; i++) {
                        assert(first + i < last);
                        operator[](i) = *(first+i);
                }
        }

        exponent_t<S, T>& operator*=(const exponent_t<S, T>& exponent) {
                for (size_t i = 0; i < S; i++) {
                        operator[](i) += exponent[i];
                }
                return *this;
        }
        exponent_t<S, T>& operator/=(const exponent_t<S, T>& exponent) {
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
        double log_eval(const boost::array<double, S>& val) const {
                double result = 0.0;
                for (size_t i = 0; i < S; i++) {
                        result += operator[](i)*log(val[i]);
                }
                return result;
        }
        using boost::array<T, S>::operator[];
};

template <size_t S, typename T, typename C = double>
class polynomial_term_t : public std::pair<exponent_t<S, T>, C> {
public:
        polynomial_term_t()
                : std::pair<exponent_t<S, T>, C>(exponent_t<S, T>(), 0.0)
                { }
        polynomial_term_t(const C& constant)
                : std::pair<exponent_t<S, T>, C>(exponent_t<S, T>(), constant)
                { }
        polynomial_term_t(const std::pair<exponent_t<S, T>, C>& pair)
                : std::pair<exponent_t<S, T>, C>(pair)
                { }
        exponent_t<S, T>& exponent() {
                return std::pair<exponent_t<S, T>, C>::first;
        }
        const exponent_t<S, T>& exponent() const {
                return std::pair<exponent_t<S, T>, C>::first;
        }
        C& coefficient() {
                return std::pair<exponent_t<S, T>, C>::second;
        }
        const C& coefficient() const {
                return std::pair<exponent_t<S, T>, C>::second;
        }

        polynomial_term_t<S, T, C> operator-() const {
                return std::pair<exponent_t<S, T>, C>(exponent(), -coefficient());
        }
        polynomial_term_t<S, T, C>& operator+=(const polynomial_term_t<S, T, C>& term) {
                coefficient() += term.coefficient();
                return *this;
        }
        polynomial_term_t<S, T, C>& operator-=(const polynomial_term_t<S, T, C>& term) {
                coefficient() -= term.coefficient();
                return *this;
        }
        polynomial_term_t<S, T, C>& operator*=(const C& constant) {
                coefficient() *= constant;
                return *this;
        }
        polynomial_term_t<S, T, C>& operator/=(const C& constant) {
                coefficient() /= constant;
                return *this;
        }
        polynomial_term_t<S, T, C>& operator*=(const polynomial_term_t<S, T, C>& term) {
                coefficient() *= term.coefficient();
                exponent()    *= term.exponent();
                return *this;
        }
        polynomial_term_t<S, T, C>& operator/=(const polynomial_term_t<S, T, C>& term) {
                coefficient() /= term.coefficient();
                exponent()    /= term.exponent();
                return *this;
        }
        C eval(const boost::array<double, S>& val) const {
                C c(coefficient());
                c *= exponent().eval(val);
                return c;
        }
        C log_eval(const boost::array<double, S>& val) const {
                C c(log(coefficient()));
                c += exponent().log_eval(val);
                return c;
        }
};

template <size_t S, typename T, typename C>
polynomial_term_t<S, T, C> operator+(const polynomial_term_t<S, T, C>& term1, const polynomial_term_t<S, T, C>& term2) {
        polynomial_term_t<S, T, C> tmp(term1);
        tmp += term2;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_term_t<S, T, C> operator-(const polynomial_term_t<S, T, C>& term1, const polynomial_term_t<S, T, C>& term2) {
        polynomial_term_t<S, T, C> tmp(term1);
        tmp -= term2;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_term_t<S, T, C> operator*(const polynomial_term_t<S, T, C>& term, const C& constant) {
        polynomial_term_t<S, T, C> tmp(term);
        tmp *= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_term_t<S, T, C> operator/(const polynomial_term_t<S, T, C>& term, const C& constant) {
        polynomial_term_t<S, T, C> tmp(term);
        tmp /= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_term_t<S, T, C> operator*(const C& constant, const polynomial_term_t<S, T, C>& term) {
        polynomial_term_t<S, T, C> tmp(term);
        tmp *= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_term_t<S, T, C> operator/(const C& constant, const polynomial_term_t<S, T, C>& term) {
        polynomial_term_t<S, T, C> tmp(constant);
        tmp /= term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_term_t<S, T, C> operator*(const polynomial_term_t<S, T, C>& term1, const polynomial_term_t<S, T, C>& term2) {
        polynomial_term_t<S, T, C> tmp(term1);
        tmp *= term2;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_term_t<S, T, C> operator/(const polynomial_term_t<S, T, C>& term1, const polynomial_term_t<S, T, C>& term2) {
        polynomial_term_t<S, T, C> tmp(term1);
        tmp /= term2;
        return tmp;
}

template <size_t S, typename T, typename C = double>
class polynomial_t : public boost::unordered_map<exponent_t<S, T>, C> {
public:
        polynomial_t()
                : boost::unordered_map<exponent_t<S, T>, C>() {}
        polynomial_t(const C& constant)
                : boost::unordered_map<exponent_t<S, T>, C>() {
                this->operator+=(polynomial_term_t<S, T, C>(constant));
        }
        polynomial_t(const exponent_t<S, T>& exponent)
                : boost::unordered_map<exponent_t<S, T>, C>() {
                this->operator+=(polynomial_term_t<S, T, C>(std::pair<exponent_t<S, T>, C>(exponent, 1.0)));
        }
        polynomial_t(const polynomial_term_t<S, T, C>& term)
                : boost::unordered_map<exponent_t<S, T>, C>() {
                this->operator+=(term);
        }

        polynomial_t<S, T, C>& operator+=(const C& constant) {
                this->operator+=(polynomial_term_t<S, T, C>(constant));
                return *this;
        }
        polynomial_t<S, T, C>& operator-=(const C& constant) {
                this->operator-=(polynomial_term_t<S, T, C>(constant));
                return *this;
        }
        polynomial_t<S, T, C>& operator+=(const polynomial_term_t<S, T, C>& term) {
                if (term.coefficient() != 0.0) {
                        C& coefficient = operator[](term.exponent());
                        coefficient += term.coefficient();
                        if (coefficient == 0.0) {
                                erase(term.exponent());
                        }
                }
                return *this;
        }
        polynomial_t<S, T, C>& operator-=(const polynomial_term_t<S, T, C>& term) {
                if (term.coefficient() != 0.0) {
                        C& coefficient = operator[](term.exponent());
                        coefficient -= term.coefficient();
                        if (coefficient == 0.0) {
                                erase(term.exponent());
                        }
                }
                return *this;
        }
        polynomial_t<S, T, C>& operator+=(const polynomial_t<S, T, C>& poly) {
                for (const_iterator it = poly.begin(); it != poly.end(); it++) {
                        operator+=(*it);
                }
                return *this;
        }
        polynomial_t<S, T, C>& operator-=(const polynomial_t<S, T, C>& poly) {
                for (const_iterator it = poly.begin(); it != poly.end(); it++) {
                        operator-=(*it);
                }
                return *this;
        }
        polynomial_t<S, T, C>& operator*=(const C& constant) {
                polynomial_t<S, T, C> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) * constant;
                }
                this->operator=(tmp);

                return *this;
        }
        polynomial_t<S, T, C>& operator/=(const C& constant) {
                polynomial_t<S, T, C> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) / constant;
                }
                this->operator=(tmp);

                return *this;
        }
        polynomial_t<S, T, C>& operator*=(const polynomial_term_t<S, T, C>& term) {
                polynomial_t<S, T, C> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) * term;
                }
                this->operator=(tmp);

                return *this;
        }
        polynomial_t<S, T, C>& operator/=(const polynomial_term_t<S, T, C>& term) {
                polynomial_t<S, T, C> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) / term;
                }
                this->operator=(tmp);

                return *this;
        }
        polynomial_t<S, T, C>& operator*=(const polynomial_t<S, T, C>& poly) {
                polynomial_t<S, T, C> tmp;

                for (typename polynomial_t<S, T, C>::const_iterator it = this->begin(); it != this->end(); it++) {
                        for (typename polynomial_t<S, T, C>::const_iterator is = poly.begin(); is != poly.end(); is++) {
                                tmp += (*it)*(*is);
                        }
                }
                this->operator=(tmp);

                return *this;
        }
        C eval(const boost::array<double, S>& val) const {
                C result(0.0);
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        result += it->eval(val);
                }
                return result;
        }
        C log_eval(const boost::array<double, S>& val) const {
                C result(-HUGE_VAL);
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        result = logadd(result, it->log_eval(val));
                }
                return result;
        }
        /* normalize coefficients so they add up to one */
        polynomial_t<S, T, C> normalize(void) const {
                polynomial_t<S, T, C> result;
                double norm = 0;

                /* compute normalization constant */
                for (typename polynomial_t<S, T, C>::const_iterator it = this->begin(); it != this->end(); it++)
                {
                        norm += it->coefficient();
                }
                for (typename polynomial_t<S, T, C>::const_iterator it = this->begin(); it != this->end(); it++)
                {
                        result += (*it) / norm;
                }
                return result;
        }

        using boost::unordered_map<exponent_t<S, T>, C>::operator[];
        using boost::unordered_map<exponent_t<S, T>, C>::erase;

        // Iterator
        ////////////////////////////////////////////////////////////////////////
        class const_iterator : public boost::unordered_map<exponent_t<S, T>, C>::const_iterator
        {
        public:
                const_iterator(typename boost::unordered_map<exponent_t<S, T>, C>::const_iterator iterator)
                        : boost::unordered_map<exponent_t<S, T>, C>::const_iterator(iterator)
                        {}

                const polynomial_term_t<S, T, C>* operator->() const
                {
                        return (const polynomial_term_t<S, T, C>*)boost::unordered_map<exponent_t<S, T>, C>::const_iterator::operator->();
                }
                const polynomial_term_t<S, T, C> operator*() const
                {
                        return *operator->();
                }
        };
        const_iterator begin() const {
                return const_iterator(boost::unordered_map<exponent_t<S, T>, C>::begin());
        }
};

template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator+(const polynomial_t<S, T, C>& poly, const C& constant) {
        polynomial_t<S, T, C> tmp(poly);
        tmp += constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator-(const polynomial_t<S, T, C>& poly, const C& constant) {
        polynomial_t<S, T, C> tmp(poly);
        tmp -= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator+(const C& constant, const polynomial_t<S, T, C>& poly) {
        polynomial_t<S, T, C> tmp(poly);
        tmp += constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator-(const C& constant, const polynomial_t<S, T, C>& poly) {
        polynomial_t<S, T, C> tmp(poly);
        tmp -= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator*(const polynomial_t<S, T, C>& poly, const C& constant) {
        polynomial_t<S, T, C> tmp(poly);
        tmp *= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator/(const polynomial_t<S, T, C>& poly, const C& constant) {
        polynomial_t<S, T, C> tmp(poly);
        tmp /= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator*(const C& constant, const polynomial_t<S, T, C>& poly) {
        polynomial_t<S, T, C> tmp(poly);
        tmp *= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator/(const C& constant, const polynomial_t<S, T, C>& poly) {
        polynomial_t<S, T, C> tmp(poly);
        tmp /= constant;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator+(const polynomial_t<S, T, C>& poly, const polynomial_term_t<S, T, C>& term) {
        polynomial_t<S, T, C> tmp(poly);
        tmp += term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator-(const polynomial_t<S, T, C>& poly, const polynomial_term_t<S, T, C>& term) {
        polynomial_t<S, T, C> tmp(poly);
        tmp -= term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator+(const polynomial_term_t<S, T, C>& term, const polynomial_t<S, T, C>& poly) {
        polynomial_t<S, T, C> tmp(poly);
        tmp += term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator-(const polynomial_term_t<S, T, C>& term, const polynomial_t<S, T, C>& poly) {
        polynomial_t<S, T, C> tmp(poly);
        tmp -= term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator*(const polynomial_t<S, T, C>& poly, const polynomial_term_t<S, T, C>& term) {
        polynomial_t<S, T, C> tmp(poly);
        tmp *= term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator/(const polynomial_t<S, T, C>& poly, const polynomial_term_t<S, T, C>& term) {
        polynomial_t<S, T, C> tmp(poly);
        tmp /= term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator*(const polynomial_term_t<S, T, C>& term, const polynomial_t<S, T, C>& poly) {
        polynomial_t<S, T, C> tmp(poly);
        tmp *= term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator/(const polynomial_term_t<S, T, C>& term, const polynomial_t<S, T, C>& poly) {
        polynomial_t<S, T, C> tmp(poly);
        tmp /= term;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator+(const polynomial_t<S, T, C>& poly1, const polynomial_t<S, T, C>& poly2) {
        polynomial_t<S, T, C> tmp(poly1);
        tmp += poly2;
        return tmp;
}
template <size_t S, typename T, typename C>
polynomial_t<S, T, C> operator*(const polynomial_t<S, T, C>& poly1, const polynomial_t<S, T, C>& poly2) {
        if (poly1.size() == 0 || poly2.size() == 0) {
                return polynomial_t<S, T, C>();
        }
        polynomial_t<S, T, C> tmp(poly1);
        tmp *= poly2;
        return tmp;
}

#endif /* _POLYNOMIAL_H_ */
