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

#ifndef __BOOST_MATH_DISTRIBUTIONS_DIRICHLET_HPP__
#define __BOOST_MATH_DISTRIBUTIONS_DIRICHLET_HPP__

#include <cmath>
#include <vector>

#include <tfbayes/utility/multinomial-beta.hh>
#include <tfbayes/utility/summation.hh>

namespace boost { namespace math {

template <class RealType = double,
          class Policy   = policies::policy<> >
class dirichlet_distribution
{
public:
        typedef RealType  input_type;
        typedef RealType result_type;

        template <class T>
        dirichlet_distribution(size_t k, T alpha)
                : m_alpha(k, alpha)
                { }
        template <class T>
        dirichlet_distribution(const T& alpha)
                : m_alpha(alpha.begin(), alpha.end())
                { }

        const std::vector<input_type>& alpha() const {
                return m_alpha;
        }
private:
        std::vector<input_type> m_alpha;
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

#endif /* __BOOST_MATH_DISTRIBUTIONS_DIRICHLET_HPP__ */
