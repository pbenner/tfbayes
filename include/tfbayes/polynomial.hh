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

template <typename T, size_t S>
class polynomial_term_t : public std::pair<exponent_t<T, S>, double> {
public:
        polynomial_term_t()
                : std::pair<exponent_t<T, S>, double>(exponent_t<T, S>(), 1.0)
                { }
        polynomial_term_t(const double constant)
                : std::pair<exponent_t<T, S>, double>(exponent_t<T, S>(), constant)
                { }
        polynomial_term_t(const std::pair<exponent_t<T, S>, double>& pair)
                : std::pair<exponent_t<T, S>, double>(pair)
                { }
        exponent_t<T, S>& exponent() {
                return std::pair<exponent_t<T, S>, double>::first;
        }
        const exponent_t<T, S>& exponent() const {
                return std::pair<exponent_t<T, S>, double>::first;
        }
        double& coefficient() {
                return std::pair<exponent_t<T, S>, double>::second;
        }
        const double& coefficient() const {
                return std::pair<exponent_t<T, S>, double>::second;
        }

        polynomial_term_t<T, S>& operator+=(const polynomial_term_t<T, S>& term) {
                coefficient() += term.coefficient();
                return *this;
        }
        polynomial_term_t<T, S>& operator-=(const polynomial_term_t<T, S>& term) {
                coefficient() -= term.coefficient();
                return *this;
        }
        polynomial_term_t<T, S>& operator*=(double constant) {
                coefficient() *= constant;
                return *this;
        }
        polynomial_term_t<T, S>& operator/=(double constant) {
                coefficient() /= constant;
                return *this;
        }
        polynomial_term_t<T, S>& operator*=(const polynomial_term_t<T, S>& term) {
                coefficient() *= term.coefficient();
                exponent()    *= term.exponent();
                return *this;
        }
        polynomial_term_t<T, S>& operator/=(const polynomial_term_t<T, S>& term) {
                coefficient() /= term.coefficient();
                exponent()    /= term.exponent();
                return *this;
        }
        double eval(const boost::array<double, S>& val) const {
                return coefficient()*exponent().eval(val);
        }
};

template <typename T, size_t S>
polynomial_term_t<T, S> operator+(const polynomial_term_t<T, S>& term1, const polynomial_term_t<T, S>& term2) {
        polynomial_term_t<T, S> tmp(term1);
        tmp += term2;
        return tmp;
}
template <typename T, size_t S>
polynomial_term_t<T, S> operator-(const polynomial_term_t<T, S>& term1, const polynomial_term_t<T, S>& term2) {
        polynomial_term_t<T, S> tmp(term1);
        tmp -= term2;
        return tmp;
}
template <typename T, size_t S>
polynomial_term_t<T, S> operator*(const polynomial_term_t<T, S>& term, double constant) {
        polynomial_term_t<T, S> tmp(term);
        tmp *= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_term_t<T, S> operator/(const polynomial_term_t<T, S>& term, double constant) {
        polynomial_term_t<T, S> tmp(term);
        tmp /= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_term_t<T, S> operator*(double constant, const polynomial_term_t<T, S>& term) {
        polynomial_term_t<T, S> tmp(term);
        tmp *= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_term_t<T, S> operator/(double constant, const polynomial_term_t<T, S>& term) {
        polynomial_term_t<T, S> tmp(term);
        tmp /= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_term_t<T, S> operator*(const polynomial_term_t<T, S>& term1, const polynomial_term_t<T, S>& term2) {
        polynomial_term_t<T, S> tmp(term1);
        tmp *= term2;
        return tmp;
}
template <typename T, size_t S>
polynomial_term_t<T, S> operator/(const polynomial_term_t<T, S>& term1, const polynomial_term_t<T, S>& term2) {
        polynomial_term_t<T, S> tmp(term1);
        tmp /= term2;
        return tmp;
}

template <typename T, size_t S>
class polynomial_t : public boost::unordered_map<exponent_t<T, S>, double> {
public:
        polynomial_t()
                : boost::unordered_map<exponent_t<T, S>, double>(),
                  _constant(0.0) {}
        polynomial_t(double constant)
                : boost::unordered_map<exponent_t<T, S>, double>(),
                  _constant(constant) {}
        polynomial_t(const polynomial_term_t<T, S>& term)
                : boost::unordered_map<exponent_t<T, S>, double>(),
                  _constant(0.0) {
                operator+=(term);
        }

        double& constant() {
                return _constant;
        }
        const double& constant() const {
                return _constant;
        }

        polynomial_t<T, S>& operator+=(double constant) {
                _constant += constant;
                return *this;
        }
        polynomial_t<T, S>& operator-=(double constant) {
                _constant -= constant;
                return *this;
        }
        polynomial_t<T, S>& operator+=(const polynomial_term_t<T, S>& term) {
                if (term.coefficient() != 0.0) {
                        operator[](term.exponent()) += term.coefficient();
                        if (operator[](term.exponent()) == 0.0) {
                                erase(term.exponent());
                        }
                }
                return *this;
        }
        polynomial_t<T, S>& operator-=(const polynomial_term_t<T, S>& term) {
                if (term.coefficient() != 0.0) {
                        operator[](term.exponent()) -= term.coefficient();
                        if (operator[](term.exponent()) == 0.0) {
                                erase(term.exponent());
                        }
                }
                return *this;
        }
        polynomial_t<T, S>& operator+=(const polynomial_t<T, S>& poly) {
                for (const_iterator it = poly.begin(); it != poly.end(); it++) {
                        operator+=(*it);
                }
                _constant += poly.constant();
                return *this;
        }
        polynomial_t<T, S>& operator-=(const polynomial_t<T, S>& poly) {
                for (const_iterator it = poly.begin(); it != poly.end(); it++) {
                        operator-=(*it);
                }
                _constant -= poly.constant();
                return *this;
        }
        polynomial_t<T, S>& operator*=(double constant) {
                polynomial_t<T, S> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) * constant;
                }
                tmp += _constant*constant;
                operator=(tmp);

                return *this;
        }
        polynomial_t<T, S>& operator/=(double constant) {
                polynomial_t<T, S> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) / constant;
                }
                tmp += _constant/constant;
                operator=(tmp);

                return *this;
        }
        polynomial_t<T, S>& operator*=(const polynomial_term_t<T, S>& term) {
                polynomial_t<T, S> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) * term;
                }
                tmp += _constant*term;
                operator=(tmp);

                return *this;
        }
        polynomial_t<T, S>& operator/=(const polynomial_term_t<T, S>& term) {
                polynomial_t<T, S> tmp;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        tmp += (*it) / term;
                }
                tmp += _constant/term;
                operator=(tmp);

                return *this;
        }
        polynomial_t<T, S>& operator*=(const polynomial_t<T, S>& poly) {
                polynomial_t<T, S> tmp;

                for (typename polynomial_t<T, S>::const_iterator it = this->begin(); it != this->end(); it++) {
                        for (typename polynomial_t<T, S>::const_iterator is = poly.begin(); is != poly.end(); is++) {
                                tmp += (*it)*(*is);
                        }
                        tmp += (*it)*poly.constant();
                }
                if (constant() != 0.0) {
                        for (typename polynomial_t<T, S>::const_iterator is = poly.begin(); is != poly.end(); is++) {
                                tmp += constant()*(*is);
                        }
                        tmp += constant()*poly.constant();
                }
                operator=(tmp);
        
                return *this;
        }
        double eval(const boost::array<double, S>& val) const {
                double result = _constant;
                for (const_iterator it = this->begin(); it != this->end(); it++) {
                        result += it->eval(val);
                }
                return result;
        }
        using boost::unordered_map<exponent_t<T, S>, double>::operator[];
        using boost::unordered_map<exponent_t<T, S>, double>::erase;

        // Iterator
        ////////////////////////////////////////////////////////////////////////
        class const_iterator : public boost::unordered_map<exponent_t<T, S>, double>::const_iterator
        {
        public:
                const_iterator(typename boost::unordered_map<exponent_t<T, S>, double>::const_iterator iterator)
                        : boost::unordered_map<exponent_t<T, S>, double>::const_iterator(iterator)
                        {}

                const polynomial_term_t<T, S>* operator->() const
                {
                        return (const polynomial_term_t<T, S>*)boost::unordered_map<exponent_t<T, S>, double>::const_iterator::operator->();
                }
                const polynomial_term_t<T, S> operator*() const
                {
                        return *operator->();
                }
        };
        const_iterator begin() const {
                return const_iterator(boost::unordered_map<exponent_t<T, S>, double>::begin());
        }

protected:
        double _constant;
};

template <typename T, size_t S>
polynomial_t<T, S> operator+(const polynomial_t<T, S>& poly, double constant) {
        polynomial_t<T, S> tmp(poly);
        tmp += constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator-(const polynomial_t<T, S>& poly, double constant) {
        polynomial_t<T, S> tmp(poly);
        tmp -= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator+(double constant, const polynomial_t<T, S>& poly) {
        polynomial_t<T, S> tmp(poly);
        tmp += constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator-(double constant, const polynomial_t<T, S>& poly) {
        polynomial_t<T, S> tmp(poly);
        tmp -= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator*(const polynomial_t<T, S>& poly, double constant) {
        polynomial_t<T, S> tmp(poly);
        tmp *= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator/(const polynomial_t<T, S>& poly, double constant) {
        polynomial_t<T, S> tmp(poly);
        tmp /= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator*(double constant, const polynomial_t<T, S>& poly) {
        polynomial_t<T, S> tmp(poly);
        tmp *= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator/(double constant, const polynomial_t<T, S>& poly) {
        polynomial_t<T, S> tmp(poly);
        tmp /= constant;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator+(const polynomial_t<T, S>& poly, const polynomial_term_t<T, S>& term) {
        polynomial_t<T, S> tmp(poly);
        tmp += term;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator-(const polynomial_t<T, S>& poly, const polynomial_term_t<T, S>& term) {
        polynomial_t<T, S> tmp(poly);
        tmp -= term;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator+(const polynomial_term_t<T, S>& term, const polynomial_t<T, S>& poly) {
        polynomial_t<T, S> tmp(poly);
        tmp += term;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator-(const polynomial_term_t<T, S>& term, const polynomial_t<T, S>& poly) {
        polynomial_t<T, S> tmp(poly);
        tmp -= term;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator*(const polynomial_t<T, S>& poly, const polynomial_term_t<T, S>& term) {
        polynomial_t<T, S> tmp(poly);
        tmp *= term;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator/(const polynomial_t<T, S>& poly, const polynomial_term_t<T, S>& term) {
        polynomial_t<T, S> tmp(poly);
        tmp /= term;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator*(const polynomial_term_t<T, S>& term, const polynomial_t<T, S>& poly) {
        polynomial_t<T, S> tmp(poly);
        tmp *= term;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator/(const polynomial_term_t<T, S>& term, const polynomial_t<T, S>& poly) {
        polynomial_t<T, S> tmp(poly);
        tmp /= term;
        return tmp;
}
template <typename T, size_t S>
polynomial_t<T, S> operator*(const polynomial_t<T, S>& poly1, const polynomial_t<T, S>& poly2) {
        polynomial_t<T, S> tmp(poly1);
        tmp *= poly2;
        return tmp;
}

#endif /* _POLYNOMIAL_H_ */
