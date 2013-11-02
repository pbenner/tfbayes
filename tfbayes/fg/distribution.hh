/* Copyright (C) 2013 Philipp Benner
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

#ifndef FG_DISTRIBUTION_HH
#define FG_DISTRIBUTION_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#define _USE_MATH_DEFINES
#include <cmath>

#include <utility/clonable.hh>

#include <boost/array.hpp>
#include <boost/noncopyable.hpp>

// interface for exponential families
class distribution_i : public clonable, boost::noncopyable {
public:
        virtual distribution_i* clone() const = 0;

        template<size_t K> double moment() const {
                return 0.0;
        }
protected:
        virtual double moment_first () const = 0;
        virtual double moment_second() const = 0;
};

class exponential_family_i : public distribution_i, boost::noncopyable {
public:

        virtual double density(double x) const = 0;
        virtual double base_measure(double x) const = 0;
        virtual double log_partition() const = 0;

        virtual exponential_family_i& operator*=(const exponential_family_i& e) = 0;
};

template<> double distribution_i::moment<1>() const {
        return moment_first();
}
template<> double distribution_i::moment<2>() const {
        return moment_second();
}

// implement what is common to most exponential families
template <size_t D>
class exponential_family_t : public exponential_family_i {
public:
        // typedefs
        typedef boost::array<double, D> array_t;

        // constructors
        exponential_family_t() { }
        exponential_family_t(const exponential_family_t& e)
        : _parameters(e.parameters()) {
        }

        friend void swap(exponential_family_t& left,
                         exponential_family_t& right) {
                using std::swap;
                swap(left._parameters, right._parameters);
        }

        // clone object
        virtual exponential_family_t* clone() const = 0;

        // pure functions
        virtual array_t statistics(double x) const = 0;

        // methods
        virtual double density(double x) const {
                array_t T = statistics(x);
                double tmp = 0.0;
                for (size_t i = 0; i < D; i++) {
                        tmp += parameters()[i]*T[i];
                }
                return base_measure(x)*std::exp(tmp - log_partition());
        }
        virtual const array_t& parameters() const {
                return _parameters;
        }
        virtual exponential_family_t& operator*=(const exponential_family_i& _e) {
                const exponential_family_t& e =
                        static_cast<const exponential_family_t&>(_e);
                for (size_t i = 0; i < D; i++) {
                        _parameters[i] += e.parameters()[i];
                }
                return *this;
        }
protected:
        virtual array_t& parameters() {
                return _parameters;
        }
        // the parameters include the base measure, which makes
        // multiplication of multiple distributions much easier
        array_t _parameters;
};

class normal_distribution_t : public exponential_family_t<2> {
public:
        normal_distribution_t(double mean, double precision)
                : _n(1.0) {
                parameters()[0] = mean*precision;
                parameters()[1] = -0.5*precision;
        }
        normal_distribution_t(const normal_distribution_t& normal_distribution)
                : exponential_family_t(normal_distribution),
                  _n      (normal_distribution._n) {
        }

        virtual normal_distribution_t* clone() const {
                return new normal_distribution_t(*this);
        }

        friend void swap(normal_distribution_t& left,
                         normal_distribution_t& right) {
                using std::swap;
                swap(static_cast<exponential_family_t<2>&>(left),
                     static_cast<exponential_family_t<2>&>(right));
                swap(left._n, right._n);
        }

        virtual double base_measure(double x) const {
                return 1.0/std::sqrt(2.0*M_PI);
        }
        virtual double log_partition() const {
                const double& p1 = parameters()[0];
                const double& p2 = parameters()[1];
                return 1.0/4.0*std::pow(p1/p2, 2.0)*p2 -
                        std::log(std::sqrt(-2*p2));
        }
        virtual array_t statistics(double x) const {
                array_t T;
                T[0] = x;
                T[1] = std::pow(x, 2.0);
                return T;
        }
protected:
        virtual double moment_first () const {
                return -0.5*parameters()[0]/parameters()[1];
        }
        virtual double moment_second() const {
                return -0.5/parameters()[1];
        }
        double _n;
};

#endif /* FG_DISTRIBUTION_HH */
