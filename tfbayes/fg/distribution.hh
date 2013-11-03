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
class distribution_i : public clonable {
public:
        virtual distribution_i* clone() const = 0;

        template<size_t K> double moment() const {
                return std::numeric_limits<double>::infinity();
        }
protected:
        virtual double moment_first () const = 0;
        virtual double moment_second() const = 0;
};

class exponential_family_i : public distribution_i {
public:

        virtual double density(double x) const = 0;
        virtual double base_measure(double x) const = 0;
        virtual double log_partition() const = 0;
        virtual void renormalize() = 0;

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
        exponential_family_t()
                : _log_partition(0.0),
                  _n            (1.0)
                { }
        exponential_family_t(const exponential_family_t& e)
                : _parameters   (e._parameters),
                  _log_partition(e._log_partition),
                  _n            (e._n)
                { }

        friend void swap(exponential_family_t& left,
                         exponential_family_t& right) {
                using std::swap;
                swap(left._parameters,    right._parameters);
                swap(left._log_partition, right._log_partition);
                swap(left._n,             right._n);
        }

        // clone object
        virtual exponential_family_t* clone() const = 0;

        // pure functions
        virtual array_t statistics(double x) const = 0;

        // methods
        virtual double density(double x) const {
                array_t T = statistics(x);
                double tmp = 0.0;
                // compute dot product
                for (size_t i = 0; i < D; i++) {
                        tmp += parameters()[i]*T[i];
                }
                return base_measure(x)*std::exp(tmp - log_partition());
        }
        virtual const array_t& parameters() const {
                return _parameters;
        }
        virtual double log_partition() const {
                return _log_partition;
        }
        virtual exponential_family_t& operator*=(const exponential_family_i& _e) {
                const exponential_family_t& e =
                        static_cast<const exponential_family_t&>(_e);
                // add parameters
                for (size_t i = 0; i < D; i++) {
                        _parameters[i] += e.parameters()[i];
                }
                // add normalization constants
                _log_partition += e._log_partition;
                // increment number of multiplications
                _n += e._n;
                // and return
                return *this;
        }
        virtual void renormalize() {
                _n = 1.0;
        }
protected:
        // access methods
        virtual array_t& parameters() {
                return _parameters;
        }
        virtual const double& n() const {
                return _n;
        }
        virtual double& n() {
                return _n;
        }
        array_t _parameters;
        // normalization constant
        double _log_partition;
        // keep track of the number of multiplications
        double _n;
};

class normal_distribution_t : public exponential_family_t<2> {
public:
        normal_distribution_t(double mean, double precision) {
                parameters()[0] = mean*precision;
                parameters()[1] = -0.5*precision;
                renormalize();
        }
        normal_distribution_t(const normal_distribution_t& normal_distribution)
                : exponential_family_t<2>(normal_distribution)
                { }

        virtual normal_distribution_t* clone() const {
                return new normal_distribution_t(*this);
        }

        friend void swap(normal_distribution_t& left,
                         normal_distribution_t& right) {
                using std::swap;
                swap(static_cast<exponential_family_t<2>&>(left),
                     static_cast<exponential_family_t<2>&>(right));
        }

        virtual double base_measure(double x) const {
                return 1.0/std::pow(std::sqrt(2.0*M_PI), n());
        }
        virtual array_t statistics(double x) const {
                array_t T;
                T[0] = x;
                T[1] = std::pow(x, 2.0);
                return T;
        }
        virtual void renormalize() {
                exponential_family_t<2>::renormalize();
                const double& p1 = parameters()[0];
                const double& p2 = parameters()[1];
                _log_partition =  -1.0/4.0*std::pow(p1/p2, 2.0)*p2 -
                        0.5*std::log(-2.0*p2);
        }
protected:
        virtual double moment_first () const {
                return -0.5*parameters()[0]/parameters()[1];
        }
        virtual double moment_second() const {
                return -0.5/parameters()[1];
        }
};

#endif /* FG_DISTRIBUTION_HH */
