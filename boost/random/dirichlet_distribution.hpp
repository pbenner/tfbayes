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

#ifndef __TFBAYES_BOOST_RANDOM_DIRICHLET_DISTRIBUTION_HPP__
#define __TFBAYES_BOOST_RANDOM_DIRICHLET_DISTRIBUTION_HPP__

#include <vector>

#include <boost/foreach.hpp>
#include <boost/random/gamma_distribution_prime.hpp>

namespace boost { namespace random {

template <class RealType = double,
          class ProbType = RealType>
class dirichlet_distribution
{
public:
        typedef RealType  input_type;
        typedef ProbType result_type;
protected:
        typedef gamma_distribution_prime<input_type, result_type> gamma_distribution_t;

        // member variables
        std::vector<gamma_distribution_t> m_rgamma;
        size_t m_size;
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
        dirichlet_distribution(const T& alpha) :
                m_size(alpha.size()) {
                BOOST_FOREACH(const typename T::value_type& a, alpha) {
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
};

} // namespace random
} // namespace boost

#endif /* __TFBAYES_BOOST_RANDOM_DIRICHLET_DISTRIBUTION_HPP__ */
