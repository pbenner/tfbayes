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

#ifndef __TFBAYES_UTILITY_DISTRIBUTION_HH__
#define __TFBAYES_UTILITY_DISTRIBUTION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/random/gamma_distribution.hpp>

#include <tfbayes/utility/probability.hh>
#include <tfbayes/utility/statistics.hh>
#include <tfbayes/utility/summation.hh>

namespace boost { namespace math {

template <class RealType, class Policy>
inline RealType pdf_derivative(const gamma_distribution<RealType, Policy>& dist, const RealType& x)
{
        RealType a = dist.shape();
        RealType b = dist.scale();

        return ((a-1)/x - 1/b)*pdf(dist, x);
}

template <class RealType, class Policy>
inline RealType log_pdf_derivative(const gamma_distribution<RealType, Policy>& dist, const RealType& x)
{
        RealType a = dist.shape();
        RealType b = dist.scale();

        return ((a-1)/x - 1/b);
}

template <class RealType = double,
          class Policy   = policies::policy<> >
class dirichlet_distribution
{
public:
        template <class T>
        dirichlet_distribution(size_t k, T alpha)
                : m_alpha(k, alpha)
                { }
        template <class T>
        dirichlet_distribution(const std::vector<T>& alpha)
                : m_alpha(alpha.begin(), alpha.end())
                { }

        const std::vector<RealType>& alpha() const {
                return m_alpha;
        }
private:
        std::vector<RealType> m_alpha;
};

template <class RealType, class Policy, class T>
inline RealType log_pdf(const dirichlet_distribution<RealType, Policy>& dist, const std::vector<T>& x)
{
        std::vector<RealType> alpha = dist.alpha();
        RealType tmp = 0.0;

        for (size_t i = 0; i < alpha.size(); i++) {
                tmp += (alpha[i]-1.0)*std::log(x[i]);
        }

        return tmp - mbeta_log(alpha);
}

template <class RealType, class Policy, class T>
inline RealType log_pdf(const dirichlet_distribution<RealType, Policy>& dist, const std::vector<T>& x, bool extended_precision)
{
        if (!extended_precision) {
                return log_pdf(dist, x);
        }
        std::vector<RealType> alpha = dist.alpha();
        std::vector<RealType> tmp(alpha.size(), 0.0);

        for (size_t i = 0; i < alpha.size(); i++) {
                tmp[i] = (alpha[i]-1.0)*std::log(x[i]);
        }

        return msum(tmp) - mbeta_log(alpha);
}

template <class RealType, class Policy, class T>
inline RealType pdf(const dirichlet_distribution<RealType, Policy>& dist, const std::vector<T>& x)
{
        return std::exp(log_pdf(dist, x));
}

template <class RealType, class Policy, class T>
inline RealType pdf(const dirichlet_distribution<RealType, Policy>& dist, const std::vector<T>& x, bool extended_precision)
{
        if (!extended_precision) {
                return log_pdf(dist, x);
        }
        return std::exp(log_pdf(dist, x));
}

} // namespace math
} // namespace boost

namespace boost { namespace random {

template <class  input_type = double,
          class result_type = input_type>
class gamma_distribution_prime
{
        input_type m_shape;
        gamma_distribution<input_type> m_gamma;
        uniform_01<input_type> m_runif;

        template<class Engine>
        result_type random_variate(Engine& eng) {
                result_type u = m_runif(eng);
                result_type b = (std::exp(1.0) + m_shape)/std::exp(1.0);
                result_type p = b*u;
                if (p <= 1.0) {
                        result_type x = std::pow(p, static_cast<result_type>(1.0/m_shape));
                        if (x == 0.0) {
                                throw std::domain_error("Invalid random variate.");
                        }
                        if (m_runif(eng) > std::exp(-x)) {
                                return random_variate(eng);
                        }
                        return x;
                }
                else {
                        result_type x = -std::log((b-p)/m_shape);
                        if (x == 0.0) {
                                throw std::domain_error("Invalid random variate.");
                        }
                        if (static_cast<result_type>(m_runif(eng)) > std::pow(x, static_cast<result_type>(m_shape-1.0))) {
                                return random_variate(eng);
                        }
                        return x;
                }
        }
public:
        gamma_distribution_prime(input_type shape)
                : m_shape(shape),
                  m_gamma(shape, 1.0)
                { }

        template<class Engine>
        result_type operator()(Engine& eng) {
                if (m_shape >= 1.0) {
                        return m_gamma(eng);
                }
                return random_variate(eng);
        }
};

template <class  input_type = double,
          class result_type = input_type>
class dirichlet_distribution
{
        typedef gamma_distribution_prime<input_type, result_type> gamma_distribution_t;
public:
        template <class T>
        dirichlet_distribution(size_t k, T alpha) :
                m_size(k) {
                for(size_t i = 0; i < k; i++) {
                        m_rgamma.push_back(
                                gamma_distribution_t(static_cast<input_type>(alpha)));
                }
        }
        template <class T>
        dirichlet_distribution(const std::vector<T>& alpha) :
                m_size(alpha.size()) {
                BOOST_FOREACH(const T& a, alpha) {
                        m_rgamma.push_back(
                                gamma_distribution_t(static_cast<input_type>(a)));
                }
        }

        template<class Engine>
        std::vector<result_type> operator()(Engine& eng) {
                result_type sum = 0.0;
                std::vector<result_type> result(m_size, 0.0);
                for (size_t i = 0; i < m_size; i++) {
                        result[i] = m_rgamma[i](eng);
                        sum += result[i];
                }
                for (size_t i = 0; i < m_size; i++) {
                        result[i] /= sum;
                }
                return result;
        }
        template<class Engine>
        std::vector<result_type> operator()(Engine& eng, bool extended_precision) {
                if (!extended_precision) {
                        return operator()(eng);
                }
                result_type sum;
                std::vector<result_type> result(m_size, 0.0);
                for (size_t i = 0; i < m_size; i++) {
                        result[i] = m_rgamma[i](eng);
                }
                sum = msum(result);
                for (size_t i = 0; i < m_size; i++) {
                        result[i] /= sum;
                }
                return result;
        }

private:
        std::vector<gamma_distribution_t> m_rgamma;
        size_t m_size;
};

} // namespace random
} // namespace boost

#endif /* __TFBAYES_UTILITY_DISTRIBUTION_HH__ */
