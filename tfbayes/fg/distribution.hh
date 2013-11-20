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

#ifndef __TFBAYES_FG_DISTRIBUTION_HH__
#define __TFBAYES_FG_DISTRIBUTION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>
#include <limits>
#include <numeric>
#define _USE_MATH_DEFINES
#include <cmath>

#include <boost/array.hpp>
#include <boost/function.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <tfbayes/fg/domain.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/default-operator.hh>
#include <tfbayes/utility/debug.hh>

// what every exponential family should provide
////////////////////////////////////////////////////////////////////////////////
class exponential_family_i : public virtual clonable {
public:
        typedef std::vector<double> vector_t;
        typedef boost::function<vector_t (const vector_t&)> statistics_f;

        virtual ~exponential_family_i() { }

        virtual exponential_family_i* clone() const = 0;

        virtual exponential_family_i& operator=(const exponential_family_i& exponential_family) = 0;

        virtual double base_measure(const vector_t& x) const = 0;
        virtual double log_partition() const = 0;
        virtual bool renormalize() = 0;
        virtual double entropy() const = 0;
        virtual vector_t moments() const = 0;
        // number of parameters
        virtual size_t k() const = 0;
        // for continuous distributions this is the density function,
        // otherwise the call operator returns a probability
        virtual double operator()(const vector_t& x) const = 0;
        virtual exponential_family_i& operator*=(const exponential_family_i& e) = 0;
        virtual bool operator==(const exponential_family_i& rhs) const = 0;
        virtual bool operator!=(const exponential_family_i& rhs) const = 0;
};
inline exponential_family_i* new_clone(const exponential_family_i& a)
{
    return a.clone();
}

// implement what is common to most exponential families
////////////////////////////////////////////////////////////////////////////////
class exponential_family_t : public exponential_family_i {
public:
        // constructors
        // n        : number of times the distribution has been multiplied
        // k        : number of parameters
        exponential_family_t(size_t k, statistics_f statistics, size_t n = 1.0,
                             const domain_t& domain = real_domain(1))
                : _parameters   (k, 0.0),
                  _statistics   (statistics),
                  _log_partition(0.0),
                  _domain       (domain),
                  _n            (n) {
        }
        exponential_family_t(const exponential_family_t& e)
                : _parameters   (e._parameters),
                  _statistics   (e._statistics),
                  _log_partition(e._log_partition),
                  _domain       (e._domain),
                  _n            (e._n)
                { }

        // clone object
        virtual exponential_family_t* clone() const = 0;

        friend void swap(exponential_family_t& left,
                         exponential_family_t& right) {
                using std::swap;
                swap(left._parameters,    right._parameters);
                swap(left._statistics,    right._statistics);
                swap(left._log_partition, right._log_partition);
                swap(left._domain,        right._domain);
                swap(left._n,             right._n);
        }

        // operators
        virtual bool operator==(const exponential_family_i& rhs) const {
                const exponential_family_t& tmp = static_cast<const exponential_family_t&>(rhs);
                for (size_t i = 0; i < k(); i++) {
                        if (std::abs(parameters()[i] - tmp.parameters()[i]) > 1.0e-10) {
                                return false;
                        }
                }
                return true;
        }
        virtual bool operator!=(const exponential_family_i& rhs) const {
                return !operator==(rhs);
        }

        // methods
        virtual double operator()(const vector_t& x) const {
                // is x in the domain of this function?
                if (!_domain.element(x)) {
                        return 0.0;
                }
                vector_t T = _statistics(x);
                // check dimensionality
                assert(T.size() == k());
                // compute density or probability
                double tmp = 0.0;
                // compute dot product
                for (size_t i = 0; i < k(); i++) {
                        tmp += parameters()[i]*T[i];
                }
                return base_measure(x)*std::exp(tmp - log_partition());
        }
        virtual const vector_t& parameters() const {
                return _parameters;
        }
        virtual double log_partition() const {
                return _log_partition;
        }
        virtual exponential_family_t& operator*=(const exponential_family_i& _e) {
                const exponential_family_t& e =
                        static_cast<const exponential_family_t&>(_e);
                // add parameters
                for (size_t i = 0; i < k(); i++) {
                        _parameters[i] += e.parameters()[i];
                }
                // add normalization constants
                _log_partition += e._log_partition;
                // increment number of multiplications
                _n += e._n;
                // recompute moments
                return *this;
        }
        virtual bool renormalize() {
                // can't renormalize if this is a null object
                if (_n == 0.0) {
                        return false;
                }
                _n = 1.0;
                return true;
        }
        virtual size_t k() const {
                return _parameters.size();
        }
protected:
        virtual vector_t& parameters() {
                return _parameters;
        }
        virtual const double& n() const {
                return _n;
        }
        virtual double& n() {
                return _n;
        }
        // natural parameters
        vector_t _parameters;
        // function computing the natural statistics
        statistics_f _statistics;
        // normalization constant
        double _log_partition;
        // domain of the density (compact support)
        domain_t _domain;
        // keep track of the number of multiplications
        double _n;
};

// a multivariate normal distribution on a simple product space
////////////////////////////////////////////////////////////////////////////////
class normal_distribution_t : public exponential_family_t {
public:
        typedef exponential_family_t base_t;

        // void object
        normal_distribution_t() :
                base_t(2, static_cast<vector_t (*)(const vector_t&)>(&statistics),
                       0.0, real_domain(1)) {
                parameters()[0] = 0.0;
                parameters()[1] = 0.0;
        }
        normal_distribution_t(double mean, double precision) :
                base_t(2, static_cast<vector_t (*)(const vector_t&)>(&statistics)
                       , 1.0, real_domain(1)) {
                assert(precision > 0.0);
                parameters()[0] = mean*precision;
                parameters()[1] = -0.5*precision;
                renormalize();
        }
        normal_distribution_t(const normal_distribution_t& normal_distribution)
                : base_t(normal_distribution)
                { }

        virtual normal_distribution_t* clone() const {
                return new normal_distribution_t(*this);
        }

        friend void swap(normal_distribution_t& left,
                         normal_distribution_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
        }
        virtual_assignment_operator(normal_distribution_t)
        derived_assignment_operator(normal_distribution_t, exponential_family_i)

        virtual double base_measure(const vector_t& x) const {
                const double m = static_cast<double>(n());
                return std::pow(2.0*M_PI, -m/2.0);
        }
        virtual bool renormalize() {
                if (!base_t::renormalize()) {
                        return false;
                }
                const double mu  = -0.5*parameters()[0]/parameters()[1];
                const double tau = -2.0*parameters()[1];
                debug("-> normal parameters:" << std::endl);
                debug("-> mean     : " << mu  << std::endl);
                debug("-> precision: " << tau << std::endl);
                _log_partition = 1.0/2.0*(tau*mu*mu - std::log(tau));
                return true;
        }
        virtual vector_t moments() const {
                vector_t m(2, 0.0);
                const double mean      = -0.5*parameters()[0]/parameters()[1];
                const double precision = -2.0*parameters()[1];
                m[0] = mean;
                m[1] = mean*mean + 1.0/precision;
                return m;
        }
        virtual double entropy() const {
                const double tau = -2.0*parameters()[1];
                return 1.0/2.0*(1.0 + std::log(2.0*M_PI*1.0/tau));
        }
        static vector_t statistics(double x);
        static vector_t statistics(const vector_t& x);
};

// the gamma distribution
////////////////////////////////////////////////////////////////////////////////
class gamma_distribution_t : public exponential_family_t {
public:
        typedef exponential_family_t base_t;

        // void object
        gamma_distribution_t() :
                base_t(2, static_cast<vector_t (*)(const vector_t&)>(&statistics),
                       0.0, positive_domain(1)) {
                parameters()[0] = 0.0;
                parameters()[1] = 0.0;
        }
        gamma_distribution_t(double shape, double rate) :
                base_t(2, static_cast<vector_t (*)(const vector_t&)>(&statistics),
                       1.0, positive_domain(1)) {
                assert(shape > 0.0);
                assert(rate  > 0.0);
                parameters()[0] =  shape-1.0;
                parameters()[1] = -rate;
                renormalize();
        }
        gamma_distribution_t(const gamma_distribution_t& gamma_distribution)
                : base_t(gamma_distribution)
                { }

        virtual gamma_distribution_t* clone() const {
                return new gamma_distribution_t(*this);
        }

        friend void swap(gamma_distribution_t& left,
                         gamma_distribution_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
        }
        virtual_assignment_operator(gamma_distribution_t)
        derived_assignment_operator(gamma_distribution_t, exponential_family_i)

        virtual double base_measure(const vector_t& x) const {
                return 1.0;
        }
        virtual bool renormalize() {
                if (!base_t::renormalize()) {
                        return false;
                }
                const double a1 =  parameters()[0]+1.0;
                const double a2 = -parameters()[1];
                debug("-> gamma parameters:" << std::endl);
                debug("-> shape: " << a1 << std::endl);
                debug("-> rate : " << a2 << std::endl);
                _log_partition =  boost::math::lgamma(a1) - (a1)*std::log(a2);
                return true;
        }
        virtual vector_t moments() const {
                vector_t m(2, 0.0);
                const double a1 =  parameters()[0]+1.0;
                const double a2 = -parameters()[1];
                m[0] = boost::math::digamma(a1) - std::log(a2);
                m[1] = a1/a2;
                return m;
        }
        virtual double entropy() const {
                const double a1 =  parameters()[0]+1.0;
                const double a2 = -parameters()[1];
                return a1 + boost::math::lgamma(a1) - std::log(a2)
                        + (1.0-a1)*boost::math::digamma(a1);
        }
        static vector_t statistics(double x);
        static vector_t statistics(const vector_t& x);
};

// the dirichlet distribution
////////////////////////////////////////////////////////////////////////////////
class dirichlet_distribution_t : public exponential_family_t {
public:
        typedef exponential_family_t base_t;

        // void object
        dirichlet_distribution_t() :
                base_t(0, &statistics, 0.0, positive_domain(1)) {
                std::fill(parameters().begin(), parameters().end(), 0.0);
        }
        dirichlet_distribution_t(const vector_t& alpha) :
                base_t(alpha.size(), &statistics, 1.0, positive_domain(alpha.size())) {
                for (size_t i = 0; i < alpha.size(); i++) {
                        assert(alpha[i] > 0.0);
                        parameters()[i] = alpha[i]-1.0;
                }
                renormalize();
        }
        dirichlet_distribution_t(const dirichlet_distribution_t& dirichlet_distribution)
                : base_t(dirichlet_distribution)
                { }

        virtual dirichlet_distribution_t* clone() const {
                return new dirichlet_distribution_t(*this);
        }

        friend void swap(dirichlet_distribution_t& left,
                         dirichlet_distribution_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
        }
        virtual_assignment_operator(dirichlet_distribution_t)
        derived_assignment_operator(dirichlet_distribution_t, exponential_family_i)

        virtual double base_measure(const vector_t& x) const {
                return 1.0;
        }
        virtual bool renormalize() {
                if (!base_t::renormalize()) {
                        return false;
                }
                double sum = std::accumulate(parameters().begin(), parameters().end(), 0.0);
                debug("-> dirichlet parameters:" << std::endl);
                _log_partition = 0.0;
                for (size_t i = 0; i < k(); i++) {
                        const double alpha = parameters()[i]+1;
                        debug(boost::format("-> alpha[%i]: %d\n") % i % alpha);
                        _log_partition += boost::math::lgamma(alpha);
                }
                _log_partition -= boost::math::lgamma(sum+k());
                return true;
        }
        virtual vector_t moments() const {
                vector_t m(k(), 0.0);
                double sum = std::accumulate(parameters().begin(), parameters().end(), 0.0);
                for (size_t i = 0; i < k(); i++) {
                        m[i] = boost::math::digamma(parameters()[i]+1.0) -
                                boost::math::digamma(sum+k());
                }
                return m;
        }
        virtual double entropy() const {
                double sum = std::accumulate(parameters().begin(), parameters().end(), 0.0);
                double h   = log_partition() + sum*boost::math::digamma(sum+k());
                for (size_t i = 0; i < k(); i++) {
                        const double alpha = parameters()[i]+1;
                        h -= (alpha-1)*boost::math::digamma(alpha);
                }
                return h;
        }
        static vector_t statistics(const vector_t& x);
};

// the discrete distribution
////////////////////////////////////////////////////////////////////////////////
class discrete_distribution_t : public exponential_family_t {
public:
        typedef exponential_family_t base_t;

        // void object
        discrete_distribution_t() :
                base_t(0, &statistics, 0.0, positive_domain(1)) {
                std::fill(parameters().begin(), parameters().end(), 0.0);
        }
        discrete_distribution_t(const vector_t& theta) :
                base_t(theta.size(), &statistics, 1.0, positive_domain(theta.size())) {
                for (size_t i = 0; i < theta.size(); i++) {
                        assert(theta[i] > 0.0);
                        parameters()[i] = std::log(theta[i]);
                }
                renormalize();
        }
        discrete_distribution_t(const discrete_distribution_t& discrete_distribution)
                : base_t(discrete_distribution)
                { }

        virtual discrete_distribution_t* clone() const {
                return new discrete_distribution_t(*this);
        }

        friend void swap(discrete_distribution_t& left,
                         discrete_distribution_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
        }
        virtual_assignment_operator(discrete_distribution_t)
        derived_assignment_operator(discrete_distribution_t, exponential_family_i)

        virtual double base_measure(const vector_t& x) const {
                return 1.0;
        }
        virtual bool renormalize() {
                if (!base_t::renormalize()) {
                        return false;
                }
                double sum = 0.0;
                debug("-> discrete parameters:" << std::endl);
                for (size_t i = 0; i < k(); i++) {
                        debug(boost::format("-> theta[%i]: %d\n") % i % std::exp(parameters()[i]));
                        sum += std::exp(parameters()[i]);
                }
                for (size_t i = 0; i < k(); i++) {
                        parameters()[i] -= std::log(sum);
                }
                _log_partition = 0.0;
                return true;
        }
        virtual vector_t moments() const {
                vector_t m(k(), 0.0);
                for (size_t i = 0; i < k(); i++) {
                        m[i] = std::exp(parameters()[i]);
                }
                return m;
        }
        virtual double entropy() const {
                return 0.0;
        }
        static vector_t statistics(const vector_t& x);
};

#endif /* __TFBAYES_FG_DISTRIBUTION_HH__ */
