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
#include <numeric>
#include <utility> /* swap */
#include <cassert>
#include <cmath>
#include <cstdint>

#include <boost/foreach.hpp>

#include <tfbayes/utility/logarithmetic.hh>

template <class RealType = double>
class probability_t
{
        typedef int8_t   sign_t;
        typedef RealType real_t;

        real_t log_p;
        sign_t sign;
public:
        probability_t()
                : log_p (-std::numeric_limits<real_t>::infinity()),
                  sign  (1)
                { }
        probability_t(real_t p)
                : log_p (std::log(std::abs(p))),
                  sign  (p >= 0.0 ? 1 : -1)
                { }
        probability_t(real_t p, sign_t sign)
                : log_p (p),
                  sign  (sign)
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
                        p  = p.abs();
                        p -= q;
                        p.sign *= -1;
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
                                p.sign *= -1;
                        }
                }
                else if (p.sign > q.sign) {
                        p += q.abs();
                }
                else {
                        p  = p.abs();
                        p += q;
                        p.sign = -1;
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
        probability_t operator-() const {
                probability_t p = *this;
                p.sign *= -1;
                return p;
        }
        bool operator==(const probability_t& p) const {
                return sign == p.sign && log_p == p.log_p;
        }
        bool operator!=(const probability_t& p) const {
                return sign != p.sign || log_p != p.log_p;
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
        explicit operator real_t() const {
                return std::exp(log_p);
        }
        explicit operator bool() const {
                return log_p != -std::numeric_limits<real_t>::infinity();
        }
        real_t log() const {
                return sign == 1
                        ? log_p
                        : std::numeric_limits<real_t>::quiet_NaN();
        }
        real_t exp() const {
                return sign == 1
                        ? std::exp( std::exp(log_p))
                        : std::exp(-std::exp(log_p));
        }
        probability_t pow(const probability_t& q) const {
                probability_t p(*this);
                assert(p.sign == 1);
                if (q.sign == 1) {
                        p.log_p *= std::exp(q.log_p);
                }
                else {
                        p.log_p /= std::exp(q.log_p);
                }
                return p;
        }
        probability_t pow(const RealType& q) const {
                probability_t p(*this);
                assert(p.sign == 1);
                p.log_p *= q;
                return p;
        }
        probability_t abs() const {
                probability_t p(*this);
                p.sign = 1;
                return p;
        }
        bool isnan() const {
                return std::isnan(log_p);
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
                if (p.log_p == -std::numeric_limits<real_t>::infinity()) {
                        o << 0.0;
                }
                else {
                        o << std::exp(p.log_p);
                }
                return o;
        }
};

template <class RealType>
probability_t<RealType> from_log_scale(RealType log_p) {
        return probability_t<RealType>(log_p, 1);
}

namespace std {
        template <class RealType>
        bool isnan(const probability_t<RealType>& p) {
                return p.isnan();
        }
        template <class RealType>
        RealType log(const probability_t<RealType>& p) {
                return p.log();
        }
        template <class RealType>
        RealType exp(const probability_t<RealType>& p) {
                return p.exp();
        }
        template <class RealType>
        probability_t<RealType> pow(const probability_t<RealType>& p, const probability_t<RealType>& q) {
                return p.pow(q);
        }
        template <class RealType>
        probability_t<RealType> pow(const probability_t<RealType>& p, const RealType& q) {
                return p.pow(q);
        }
        template <class RealType>
        probability_t<RealType> abs(const probability_t<RealType>& p) {
                return p.abs();
        }
}

template <class RealType>
RealType GCC_ATTRIBUTE_NOAMATH
kahan_sum(const std::vector<RealType>& input)
{
        RealType sum = 0.0;
        RealType c   = 0.0;
        for (size_t i = 0; i < input.size(); i++) {
                RealType y = input[i] - c;
                RealType t = sum + y;
                c = (t - sum) - y;
                sum = t;
        }
        return sum;
}

template <class RealType>
RealType GCC_ATTRIBUTE_NOAMATH
msum(const std::vector<RealType>& input) {
        // Code from: http://code.activestate.com/recipes/393090/
        // Full precision summation using multiple floats for intermediate values
        // Rounded x+y stored in hi with the round-off stored in lo.  Together
        // hi+lo are exactly equal to x+y.  The inner loop applies hi/lo summation
        // to each partial so that the list of partial sums remains exact.
        // Depends on IEEE-754 arithmetic guarantees.  See proof of correctness at:
        // www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
        RealType hi;
        RealType lo;
        std::vector<RealType> partials;
        std::vector<RealType> partials_tmp;

        BOOST_FOREACH(RealType x, input) {
                size_t i = 0;
                BOOST_FOREACH(RealType y, partials) {
                        if (std::abs(x) < std::abs(y)) {
                                std::swap(x,y);
                        }
                        hi = x + y;
                        lo = y - (hi - x);
                        if (lo) {
                                if (i < partials_tmp.size()) {
                                        partials_tmp[i] = lo;
                                }
                                else {
                                        partials_tmp.push_back(lo);
                                        i++;
                                }
                        }
                        x = hi;
                }
                while (i+1 < partials_tmp.size()) {
                        partials_tmp.pop_back();
                }
                partials_tmp.push_back(x);
                partials = partials_tmp;
        }
        return std::accumulate(partials.begin(), partials.end(), RealType(0));
}

#endif /* __TFBAYES_UTILITY_PROBABILITY_HH__ */
