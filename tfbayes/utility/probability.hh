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
#include <cmath>

#include <tfbayes/utility/logarithmetic.hh>

class probability_t;

namespace std {
        double log(const probability_t& p);
}

class probability_t {
public:
        probability_t()
        : log_p(-std::numeric_limits<double>::infinity())
                { }
        probability_t(double p)
        : log_p(std::log(p))
                { }
        probability_t(const probability_t& p)
        : log_p(p.log_p)
                { }

        probability_t& operator+=(const probability_t& p) {
                log_p = logadd(log_p, p.log_p);
                return *this;
        }
        probability_t& operator-=(const probability_t& p) {
                log_p = logsub(log_p, p.log_p);
                return *this;
        }
        probability_t& operator*=(const probability_t& p) {
                log_p += p.log_p;
                return *this;
        }
        probability_t& operator/=(const probability_t& p) {
                log_p -= p.log_p;
                return *this;
        }
        probability_t& operator%=(const probability_t& p) {
                while (log_p >= p.log_p) {
                        operator-=(p);
                }
                return *this;
        }
        bool operator>=(const probability_t& p) {
                return log_p >= p.log_p;
        }
        bool operator>(const probability_t& p) {
                return log_p > p.log_p;
        }
        bool operator<=(const probability_t& p) {
                return log_p <= p.log_p;
        }
        bool operator<(const probability_t& p) {
                return log_p < p.log_p;
        }
        explicit operator double() const {
                return std::exp(log_p);
        }

        friend
        std::ostream& operator<< (std::ostream& o, const probability_t& p) {
                if (p.log_p == -std::numeric_limits<double>::infinity()) {
                        o << 0.0;
                }
                else {
                        o << std::exp(p.log_p);
                }
                return o;
        }
        friend double std::log(const probability_t& p);

protected:
        double log_p;
};

probability_t operator+(const probability_t& p, const probability_t& q) {
        probability_t tmp(p);
        return tmp += q;
}
probability_t operator-(const probability_t& p, const probability_t& q) {
        probability_t tmp(p);
        return tmp -= q;
}
probability_t operator*(const probability_t& p, const probability_t& q) {
        probability_t tmp(p);
        return tmp *= q;
}
probability_t operator/(const probability_t& p, const probability_t& q) {
        probability_t tmp(p);
        return tmp /= q;
}
probability_t operator%(const probability_t& p, const probability_t& q) {
        probability_t tmp(p);
        return tmp %= q;
}

namespace std {

        double log(const probability_t& p) {
                return p.log_p;
        }
}

#endif /* __TFBAYES_UTILITY_PROBABILITY_HH__ */
