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
#define _USE_MATH_DEFINES
#include <cmath>

#include <boost/array.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <tfbayes/fg/domain.hh>
#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/default-operator.hh>
#include <tfbayes/utility/debug.hh>

// interface for distributions, in general a distribution might not
// even have a density (e.g. dirac measure), so this is very limited
////////////////////////////////////////////////////////////////////////////////
class distribution_i : public virtual clonable {
public:
        virtual distribution_i* clone() const = 0;

        virtual distribution_i& operator=(const distribution_i& distribution) = 0;

        // comparison operators
        virtual bool operator==(const distribution_i& rhs) const = 0;
        virtual bool operator!=(const distribution_i& rhs) const = 0;

        // dimension of the space this distribution is defined on
        virtual size_t dimension() const = 0;
};

// the simplest distribution is a dirac measure
////////////////////////////////////////////////////////////////////////////////
class dirac_distribution_t : public distribution_i {
public:
        typedef std::vector<double> x_t;

        dirac_distribution_t() :
                _x(){ }
        // one-dimensional dirac
        dirac_distribution_t(double x) :
                _x(1, x) { }
        // n-dimensional dirac
        dirac_distribution_t(const std::vector<double>& x) :
                _x(x) { }

        virtual dirac_distribution_t* clone() const {
                return new dirac_distribution_t(*this);
        }

        friend void swap(dirac_distribution_t& left,
                         dirac_distribution_t& right) {
                using std::swap;
                swap(left._x, right._x);
        }
        derived_assignment_operator(distribution_i, dirac_distribution_t)

        virtual bool operator==(const distribution_i& rhs) const {
                const dirac_distribution_t& tmp =
                        static_cast<const dirac_distribution_t&>(rhs);
                return std::equal(_x.begin(), _x.end(), tmp._x.begin());
        }
        virtual bool operator!=(const distribution_i& rhs) const {
                return !operator==(rhs);
        }
        virtual size_t dimension() const {
                return _x.size();
        }
        virtual double operator[](size_t i) const {
                return _x[i];
        }

protected:
        // location of the dirac measure
        std::vector<double> _x;
};

// a class that carries the first moments of the sufficient statistics
////////////////////////////////////////////////////////////////////////////////

class sufficient_moments_i : public virtual clonable {
public:
        virtual ~sufficient_moments_i() { }

        virtual sufficient_moments_i* clone() const = 0;

        virtual sufficient_moments_i& operator=(const sufficient_moments_i& sufficient_moments) = 0;

        virtual const double& operator[](size_t i) const = 0;
        virtual double& operator[](size_t i) = 0;

        virtual bool operator==(const sufficient_moments_i& rhs) const = 0;
        virtual bool operator!=(const sufficient_moments_i& rhs) const = 0;

        virtual size_t n() const = 0;
};

template <size_t D>
class sufficient_moments_t : public sufficient_moments_i, public boost::array<double, D> {
public:
        typedef boost::array<double, D> base_t;

        virtual sufficient_moments_t* clone() const = 0;

        virtual sufficient_moments_t& operator=(const sufficient_moments_i& sufficient_moments) = 0;

        friend void swap(sufficient_moments_t& left, sufficient_moments_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left), static_cast<base_t&>(right));
        }

        virtual const double& operator[](size_t i) const {
                return base_t::operator[](i);
        }
        virtual double& operator[](size_t i) {
                return base_t::operator[](i);
        }
        virtual bool operator==(const sufficient_moments_i& rhs) const {
                bool result = true;
                for (size_t i = 0; i < D; i++) {
                        result &= std::abs(base_t::operator[](i) - rhs[i]) < 1e-20;
                }
                return result;
        }
        virtual bool operator!=(const sufficient_moments_i& rhs) const {
                return !operator==(rhs);
        }
        virtual size_t n() const {
                return D;
        }
};

class normal_moments_t : public sufficient_moments_t<2> {
public:
        typedef sufficient_moments_t<2> base_t;
        typedef boost::array<double, 2> parameters_t;

        virtual normal_moments_t* clone() const {
                return new normal_moments_t(*this);
        }

        friend void swap(normal_moments_t& left, normal_moments_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left), static_cast<base_t&>(right));
        }
        derived_assignment_operator(sufficient_moments_i, normal_moments_t)

        normal_moments_t() { }
        normal_moments_t(double mean, double precision) {
                this->operator[](0) = mean;
                this->operator[](1) = mean*mean + 1.0/(precision*precision);
        }
        normal_moments_t(const parameters_t& parameters) {
                const double mean      = -0.5*parameters[0]/parameters[1];
                const double precision = -2.0*parameters[1];
                this->operator[](0) = mean;
                this->operator[](1) = mean*mean + 1.0/(precision*precision);
        }
        normal_moments_t(const dirac_distribution_t& dirac) {
                this->operator[](0) = 0.0;
                this->operator[](1) = 0.0;
                for (size_t i = 0; i < dirac.dimension(); i++) {
                        this->operator[](0) += dirac[i];
                        this->operator[](1) += dirac[i]*dirac[i];
                }
        }
};

class gamma_moments_t : public sufficient_moments_t<2> {
public:
        typedef sufficient_moments_t<2> base_t;
        typedef boost::array<double, 2> parameters_t;

        virtual gamma_moments_t* clone() const {
                return new gamma_moments_t(*this);
        }

        friend void swap(gamma_moments_t& left, gamma_moments_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left), static_cast<base_t&>(right));
        }
        derived_assignment_operator(sufficient_moments_i, gamma_moments_t)

        gamma_moments_t() { }
        gamma_moments_t(double shape, double rate) {
                this->operator[](0) = shape-1.0;
                this->operator[](1) = -rate;
        }
        gamma_moments_t(const parameters_t& parameters) {
                const double a1 =  parameters[0]+1.0;
                const double a2 = -parameters[1];
                this->operator[](0) = boost::math::digamma(a1) - std::log(a2);
                this->operator[](1) = a1/a2;
        }
        gamma_moments_t(const dirac_distribution_t& dirac) {
                this->operator[](0) = 0.0;
                this->operator[](1) = 0.0;
                for (size_t i = 0; i < dirac.dimension(); i++) {
                        this->operator[](0) += std::log(dirac[i]);
                        this->operator[](1) += dirac[i];
                }
        }
};

// what every exponential family should provide
////////////////////////////////////////////////////////////////////////////////
class exponential_family_i : public distribution_i {
public:
        virtual exponential_family_i* clone() const = 0;

        virtual double density(const std::vector<double>& x) const = 0;
        virtual double base_measure(const std::vector<double>& x) const = 0;
        virtual double log_partition() const = 0;
        virtual bool renormalize() = 0;

        // multiplication of exponential families
        virtual exponential_family_i& operator*=(const exponential_family_i& e) = 0;

        virtual const sufficient_moments_i& moments() const = 0;
protected:
        virtual void update_moments() = 0;
};

// implement what is common to most exponential families
////////////////////////////////////////////////////////////////////////////////
template <size_t D, typename M>
class exponential_family_t : public exponential_family_i {
public:
        // typedefs
        typedef boost::array<double, D> array_t;
        typedef M moments_t;

        // constructors
        exponential_family_t(size_t dim = 1, size_t n = 1.0, const domain_t& domain = real_domain(1))
                : _log_partition(0.0),
                  _domain       (domain),
                  _n            (n),
                  _dimension    (dim) {
                assert(dim > 0);
        }
        exponential_family_t(const exponential_family_t& e)
                : _parameters   (e._parameters),
                  _log_partition(e._log_partition),
                  _domain       (e._domain),
                  _n            (e._n),
                  _dimension    (e._dimension),
                  _moments      (e._moments)
                { }

        // clone object
        virtual exponential_family_t* clone() const = 0;

        friend void swap(exponential_family_t& left,
                         exponential_family_t& right) {
                using std::swap;
                swap(left._parameters,    right._parameters);
                swap(left._log_partition, right._log_partition);
                swap(left._domain,        right._domain);
                swap(left._n,             right._n);
                swap(left._dimension,     right._dimension);
                swap(left._moments,       right._moments);
        }

        // operators
        virtual bool operator==(const distribution_i& rhs) const {
                const exponential_family_t& tmp = static_cast<const exponential_family_t&>(rhs);
                if (dimension() != tmp.dimension()) {
                        return false;
                }
                for (size_t i = 0; i < D; i++) {
                        if (std::abs(parameters()[i] - tmp.parameters()[i]) > 1.0e-10) {
                                return false;
                        }
                }
                return true;
        }
        virtual bool operator!=(const distribution_i& rhs) const {
                return !operator==(rhs);
        }

        // pure functions
        virtual const array_t& statistics(const std::vector<double>& x) const = 0;

        // methods
        virtual double density(const std::vector<double>& x) const {
                const array_t& T = statistics(x);
                // is x in the domain of this function?
                if (!_domain.element(x)) {
                        return 0.0;
                }
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
                // recompute moments
                update_moments();
                // and return
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
        virtual size_t dimension() const {
                return _dimension;
        }
        const moments_t& moments() const {
                return _moments;
        }
protected:
        virtual void update_moments() {
                _moments = moments_t(parameters());
        }
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
        // natural parameters
        array_t _parameters;
        // sufficient statistics
        mutable array_t _T;
        // normalization constant
        double _log_partition;
        // domain of the density (compact support)
        domain_t _domain;
        // keep track of the number of multiplications
        double _n;
        // dimension of the space where this distribution is defined
        // on
        size_t _dimension;
        // moments of the sufficient statistics
        moments_t _moments;
};

// a multivariate normal distribution on a simple product space
////////////////////////////////////////////////////////////////////////////////
class normal_distribution_t : public exponential_family_t<2, normal_moments_t> {
public:
        typedef exponential_family_t<2, normal_moments_t> base_t;

        // void object
        normal_distribution_t() :
                base_t(1, 0.0, real_domain(1)) {
                parameters()[0] = 0.0;
                parameters()[1] = 0.0;
        }
        normal_distribution_t(double mean, double precision, size_t dim = 1) :
                base_t(dim, 1.0, real_domain(dim)) {
                assert(precision > 0.0);
                parameters()[0] = mean*precision;
                parameters()[1] = -0.5*precision;
                renormalize();
                update_moments();
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
        derived_assignment_operator(distribution_i, normal_distribution_t)

        virtual double base_measure(const std::vector<double>& x) const {
                return std::pow(2.0*M_PI, -static_cast<double>(dimension())/2.0*static_cast<double>(n()));
        }
        virtual const array_t& statistics(const std::vector<double>& x) const {
                assert(x.size() == dimension());
                // reset statistics
                _T[0] = 0.0;
                _T[1] = 0.0;
                for (size_t i = 0; i < dimension(); i++) {
                        _T[0] += x[i];
                        _T[1] += std::pow(x[i], 2.0);
                }
                return _T;
        }
        virtual bool renormalize() {
                if (!base_t::renormalize()) {
                        return false;
                }
                const double& mu  = -0.5*parameters()[0]/parameters()[1];
                const double& tau = -2.0*parameters()[1];
                debug("-> normal parameters:" << std::endl);
                debug("-> mean     : " << mu  << std::endl);
                debug("-> precision: " << tau << std::endl);
                _log_partition = static_cast<double>(dimension())/2.0*tau*mu*mu -
                        static_cast<double>(dimension())/2.0*std::log(tau);
                return true;
        }
};

// the gamma distribution
////////////////////////////////////////////////////////////////////////////////
class gamma_distribution_t : public exponential_family_t<2, gamma_moments_t> {
public:
        typedef exponential_family_t<2, gamma_moments_t> base_t;

        // void object
        gamma_distribution_t() :
                base_t(1, 0.0, positive_domain(1)) {
                parameters()[0] = 0.0;
                parameters()[1] = 0.0;
        }
        gamma_distribution_t(double shape, double rate, size_t dim = 1) :
                base_t(dim, 1.0, positive_domain(1)) {
                assert(shape > 0.0);
                assert(rate  > 0.0);
                parameters()[0] =  shape-1.0;
                parameters()[1] = -rate;
                renormalize();
                update_moments();
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
        derived_assignment_operator(distribution_i, gamma_distribution_t)

        virtual double base_measure(const std::vector<double>& x) const {
                return 1.0;
        }
        virtual const array_t& statistics(const std::vector<double>& x) const {
                assert(x.size() == dimension());
                // reset statistics
                _T[0] = 0.0;
                _T[1] = 0.0;
                for (size_t i = 0; i < dimension(); i++) {
                        _T[0] += std::log(x[i]);
                        _T[1] += x[i];
                }
                return _T;
        }
        virtual bool renormalize() {
                if (!base_t::renormalize()) {
                        return false;
                }
                const double& a1 = parameters()[0]+1.0;
                const double& a2 = -parameters()[1];
                debug("-> gamma parameters:" << std::endl);
                debug("-> a1: " << a1 << std::endl);
                debug("-> a2: " << a2 << std::endl);
                _log_partition =  boost::math::lgamma(a1) -
                        (a1)*std::log(a2);
                return true;
        }
};

#endif /* __TFBAYES_FG_DISTRIBUTION_HH__ */
