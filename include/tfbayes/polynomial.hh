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

#include <boost/array.hpp>
#include <boost/unordered_set.hpp>

template <typename T, size_t S>
class polynomial_term_t : public boost::array<T, S> {
public:
        polynomial_term_t()
                : _coefficient(1.0) {
                for (size_t i = 0; i < S; i++) {
                        operator[](i) = 0;
                }
        }
        double& coefficient() {
                return _coefficient;
        }
        const double coefficient() const {
                return _coefficient;
        }
        using boost::array<T, S>::operator[];

        polynomial_term_t<T, S>& operator+=(const polynomial_term_t<T, S>& term) {
                _coefficient += term.coefficient();
                return *this;
        }
        polynomial_term_t<T, S>& operator-=(const polynomial_term_t<T, S>& term) {
                _coefficient -= term.coefficient();
                return *this;
        }
        polynomial_term_t<T, S>& operator*=(double constant) {
                _coefficient *= constant;
                return *this;
        }
        polynomial_term_t<T, S>& operator/=(double constant) {
                _coefficient /= constant;
                return *this;
        }
        polynomial_term_t<T, S>& operator*=(const polynomial_term_t<T, S>& term) {
                _coefficient *= term.coefficient();
                for (size_t i = 0; i < S; i++) {
                        operator[](i) += term[i];
                }
                return *this;
        }
        polynomial_term_t<T, S>& operator/=(const polynomial_term_t<T, S>& term) {
                _coefficient /= term.coefficient();
                for (size_t i = 0; i < S; i++) {
                        operator[](i) -= term[i];
                }
                return *this;
        }

protected:
        double _coefficient;
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
class polynomial_t : public boost::unordered_set<polynomial_term_t<T, S> > {
public:
        typedef typename polynomial_t<T, S>::iterator iterator;
        typedef typename polynomial_t<T, S>::const_iterator const_iterator;

        polynomial_t()
                : boost::unordered_set<polynomial_term_t<T, S> >(),
                  _constant(0.0) {}

        double& constant() {
                return _constant;
        }
        const double constant() const {
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
                if (count(term) == 0 && term.coefficient() != 0) {
                        insert(term);
                }
                else {
                        polynomial_term_t<T, S> tmp(*find(term));
                        erase(term);
                        tmp += term;
                        if (tmp.coefficient() != 0.0) {
                                insert(tmp);
                        }
                }
                return *this;
        }
        polynomial_t<T, S>& operator-=(const polynomial_term_t<T, S>& term) {
                if (count(term) == 0 && term.coefficient() != 0) {
                        polynomial_term_t<T, S> tmp(term);
                        tmp.coefficient() = -tmp.coefficient();
                        insert(tmp);
                }
                else {
                        polynomial_term_t<T, S> tmp(*find(term));
                        erase(term);
                        tmp -= term;
                        if (tmp.coefficient() != 0.0) {
                                insert(tmp);
                        }
                }
                return *this;
        }
        polynomial_t<T, S>& operator*=(double constant) {
                for (iterator it = this->begin(); it != this->end(); it++) {
                        (*it) *= constant;
                }
                operator+=(_constant*constant);

                return *this;
        }
        polynomial_t<T, S>& operator/=(double constant) {
                for (iterator it = this->begin(); it != this->end(); it++) {
                        (*it) /= constant;
                }
                operator+=(_constant/constant);

                return *this;
        }
        polynomial_t<T, S>& operator*=(const polynomial_term_t<T, S>& term) {
                for (iterator it = this->begin(); it != this->end(); it++) {
                        (*it) *= term;
                }
                operator+=(_constant*term);

                return *this;
        }
        polynomial_t<T, S>& operator/=(const polynomial_term_t<T, S>& term) {
                for (iterator it = this->begin(); it != this->end(); it++) {
                        (*it) /= term;
                }
                operator+=(_constant/term);

                return *this;
        }
        using boost::unordered_set<polynomial_term_t<T, S> >::insert;
        using boost::unordered_set<polynomial_term_t<T, S> >::count;
        using boost::unordered_set<polynomial_term_t<T, S> >::find;

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
        polynomial_t<T, S> result;

        for (typename polynomial_t<T, S>::const_iterator it = poly1.begin(); it != poly1.end(); it++) {
                for (typename polynomial_t<T, S>::const_iterator is = poly2.begin(); is != poly2.end(); is++) {
                        result += (*it)*(*is);
                }
                result += (*it)*poly2.constant();
        }
        if (poly1.constant() != 0.0) {
                for (typename polynomial_t<T, S>::const_iterator is = poly2.begin(); is != poly2.end(); is++) {
                        result += poly1.constant()*(*is);
                }
        }
        
        return result;
}

#endif /* _POLYNOMIAL_H_ */
