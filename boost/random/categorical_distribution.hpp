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

#ifndef __TFBAYES_BOOST_RANDOM_CATEGORICAL_DISTRIBUTION_HPP__
#define __TFBAYES_BOOST_RANDOM_CATEGORICAL_DISTRIBUTION_HPP__

#include <vector>

#include <boost/foreach.hpp>
#include <boost/random/uniform_01.hpp>

namespace boost { namespace random {

template <class RealType = double,
          class  IntType = size_t>
class categorical_distribution
{
public:
        typedef RealType  input_type;
        typedef  IntType result_type;
protected:
        // vector of probabilities
        std::vector<input_type> m_p;
        // uniform distribution
        boost::random::uniform_01<input_type> runif;
public:
        categorical_distribution()
                { }
        template <class T>
        categorical_distribution(const T& alpha)
        : m_p(alpha.begin(), alpha.end())
                { }

        template<class Engine>
        result_type operator()(Engine& eng) {
        
                input_type r = runif(eng);
                input_type c = 0.0;

                for (result_type i = 0; i < m_p.size()-1; i += 1) {
                        /* comute cumulative probabilities */
                        c += m_p[i];
                        /* return if c is large enough */
                        if (r <= c) return i;
                }
                return m_p.size()-1;
        }
};

} // namespace random
} // namespace boost

#endif /* __TFBAYES_BOOST_RANDOM_CATEGORICAL_DISTRIBUTION_HPP__ */
