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

#ifndef __TFBAYES_BOOST_RANDOM_GAMMA_DISTRIBUTION_PRIME_HPP__
#define __TFBAYES_BOOST_RANDOM_GAMMA_DISTRIBUTION_PRIME_HPP__

#include <cmath>
#include <vector>

#include <boost/foreach.hpp>

#include <boost/random/gamma_distribution.hpp>

namespace boost { namespace random {

template <class RealType = double,
          class ProbType = RealType>
class gamma_distribution_prime
{
public:
        typedef RealType  input_type;
        typedef ProbType result_type;
protected:
        // member variables
        input_type m_shape;
        gamma_distribution<input_type> m_gamma;
        uniform_01<input_type> m_runif;
public:
        template<class Engine>
        result_type random_variate(Engine& eng) {
                result_type u = m_runif(eng);
                result_type b = (input_type(std::exp(1.0)) + input_type(m_shape))/input_type(std::exp(1.0));
                result_type p = b*u;
                if (p <= 1.0) {
                        result_type x = std::pow(p, 1.0/m_shape);
                        if (x == 0.0) {
                                throw std::domain_error("Invalid random variate.");
                        }
                        if (result_type(-std::log(m_runif(eng))) < x) {
                                return random_variate(eng);
                        }
                        return x;
                }
                else {
                        result_type x = -std::log((b-p)/m_shape);
                        if (x == 0.0) {
                                throw std::domain_error("Invalid random variate.");
                        }
                        if (result_type(m_runif(eng)) > std::pow(x, m_shape-1.0)) {
                                return random_variate(eng);
                        }
                        return x;
                }
        }

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


} // namespace random
} // namespace boost

#endif /* __TFBAYES_BOOST_RANDOM_GAMMA_DISTRIBUTION_PRIME_HPP__ */
