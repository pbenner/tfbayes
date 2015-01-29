/* Copyright (C) 2015 Philipp Benner
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

#ifndef __TFBAYES_UTILITY_PROBABILITY_HH__
#define __TFBAYES_UTILITY_PROBABILITY_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <utility> /* swap */
#include <cassert>
#include <cmath>
#include <cstdint>

#include <tfbayes/utility/logarithmetic.hh>

class probability_t {
public:
        probability_t()
        : log_p(-std::numeric_limits<double>::infinity())
                { }
        probability_t(double p)
                : log_p(std::log(std::abs(p))), sign(p >= 0.0 ? 1 : -1)
                { }

        probability_t& operator+=(const probability_t& q_) {
                probability_t& p = *this;
                probability_t  q = q_;
                if (p.sign == q.sign) {
                        p.log_p = logadd(p.log_p, q.log_p);
                }
                else if (p.sign > q.sign) {
                        p -= q.abs();
                }
                else {
                        p.abs(); p -= q; sign *= -1;
                }
                return *this;
        }
        probability_t& operator-=(const probability_t& q_) {
                probability_t& p = *this;
                probability_t  q = q_;
                if (p.sign == q.sign) {
                        if (p.log_p > q.log_p) {
                                p.log_p = logsub(p.log_p, q.log_p);
                        }
                        else {
                                p.log_p = logsub(q.log_p, p.log_p);
                                sign *= -1;
                        }
                }
                else if (p.sign > q.sign) {
                        p += q.abs();
                }
                else {
                        p.abs(); p += q; sign = -1;
                }
                return *this;
        }
        probability_t& operator*=(const probability_t& p) {
                log_p += p.log_p;
                sign  *= p.sign;
                return *this;
        }
        probability_t& operator/=(const probability_t& p) {
                log_p -= p.log_p;
                sign  *= p.sign;
                return *this;
        }
        // This operator behaves like the modulo operator in GNU-R!
        probability_t& operator%=(const probability_t& q_) {
                probability_t& p = *this;
                probability_t  q = q_;
                assert(q >= 0.0);
                if (p >= 0) {
                        while (p >= q) {
                                p -= q;
                        }
                } else {
                        while (p < 0.0) {
                                p += q;
                        }
                }
                return *this;
        }
        probability_t& operator-() {
                sign *= -1;
                return *this;
        }
        bool operator>=(const probability_t& p) const {
                return (sign == p.sign && log_p >= p.log_p) || sign > p.sign;
        }
        bool operator>(const probability_t& p) const {
                return (sign == p.sign && log_p > p.log_p) || sign > p.sign;
        }
        bool operator<=(const probability_t& p) const {
                return (sign == p.sign && log_p <= p.log_p) || sign < p.sign;
        }
        bool operator<(const probability_t& p) const {
                return (sign == p.sign && log_p < p.log_p) || sign < p.sign;
        }
        explicit operator double() const {
                return std::exp(log_p);
        }
        double log() const {
                return sign == 1
                        ? log_p
                        : std::numeric_limits<double>::quiet_NaN();
        }
        probability_t& abs() {
                sign = 1;
                return *this;
        }
        // define operators here so that they can be used by member functions
        friend
        probability_t operator+(const probability_t& p, const probability_t& q) {
                probability_t tmp(p);
                return tmp += q;
        }
        friend
        probability_t operator-(const probability_t& p, const probability_t& q) {
                probability_t tmp(p);
                return tmp -= q;
        }
        friend
        probability_t operator*(const probability_t& p, const probability_t& q) {
                probability_t tmp(p);
                return tmp *= q;
        }
        friend
        probability_t operator/(const probability_t& p, const probability_t& q) {
                probability_t tmp(p);
                return tmp /= q;
        }
        friend
        probability_t operator%(const probability_t& p, const probability_t& q) {
                probability_t tmp(p);
                return tmp %= q;
        }
        friend
        std::ostream& operator<< (std::ostream& o, const probability_t& p) {
                if (p.sign == -1) {
                        o << "-";
                }
                if (p.log_p == -std::numeric_limits<double>::infinity()) {
                        o << 0.0;
                }
                else {
                        o << std::exp(p.log_p);
                }
                return o;
        }

protected:
        double log_p;
        int8_t sign;
};

namespace std {
        double log(const probability_t& p) {
                return p.log();
        }
        probability_t abs(const probability_t& p) {
                probability_t tmp = p;
                return tmp.abs();
        }
}

#endif /* __TFBAYES_UTILITY_PROBABILITY_HH__ */
