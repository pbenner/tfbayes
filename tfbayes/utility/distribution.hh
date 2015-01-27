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

#include <boost/foreach.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/random/gamma_distribution.hpp>

#include <tfbayes/utility/probability.hh>
#include <tfbayes/utility/statistics.hh>

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
        dirichlet_distribution(const std::vector<RealType>& alpha)
                : m_alpha(alpha)
                { }

        const std::vector<RealType>& alpha() const {
                return m_alpha;
        }
private:
        std::vector<double> m_alpha;
};

template <class RealType, class Policy>
inline RealType log_pdf(const dirichlet_distribution<RealType, Policy>& dist, const std::vector<probability_t>& x)
{
        std::vector<RealType> alpha = dist.alpha();
        RealType tmp = 0.0;

        for (size_t i = 0; i < alpha.size(); i++) {
                tmp += (alpha[i]-1.0)*std::log(x[i]);
        }

        return tmp - mbeta_log(alpha);
}

template <class RealType, class Policy>
inline RealType pdf(const dirichlet_distribution<RealType, Policy>& dist, const std::vector<probability_t> x)
{
        return std::exp(log_pdf(dist, x));
}

} // namespace math
} // namespace boost

namespace boost { namespace random {

template <class RealType = double>
class dirichlet_distribution
{
public:
        typedef double input_type;
        typedef RealType result_type;


        dirichlet_distribution(const std::vector<input_type>& alpha) :
                m_size(alpha.size()) {
                BOOST_FOREACH(const input_type& a, alpha) {
                        m_distributions.push_back(
                                gamma_distribution<input_type>(a));
                }
        }

        template<class Engine>
        std::vector<result_type> operator()(Engine& eng) {
                RealType sum;
                std::vector<result_type> result(m_size, 0.0);
                for (size_t i = 0; i < m_size; i++) {
                        result[i] = m_distributions[i](eng);
                        sum += result[i];
                }
                for (size_t i = 0; i < m_size; i++) {
                        result[i] /= sum;
                }
                return result;
        }

private:
        std::vector<gamma_distribution<input_type> > m_distributions;
        size_t m_size;
};

} // namespace random
} // namespace boost

#endif /* __TFBAYES_UTILITY_DISTRIBUTION_HH__ */
