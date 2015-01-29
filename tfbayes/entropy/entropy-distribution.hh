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

#ifndef __TFBAYES_ENTROPY_SAMPLER_HH__
#define __TFBAYES_ENTROPY_SAMPLER_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/math/distributions/beta.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

#include <tfbayes/utility/histogram.hh>

namespace boost { namespace random {

template <class input_type = double, class result_type = input_type>
class entropy_distribution
{
public:
        entropy_distribution(size_t k, input_type a1, input_type a2) :
                m_k(k), m_a1(a1), m_a2(a2),
                m_state (k, 0.0),
                m_beta  (a1, a2) {
        }

        template<class Engine>
        const std::vector<result_type>& operator()(Engine& eng, input_type sigma = 0.01) {

                return m_state;
        }
        

private:
        result_type sum_proposal(size_t except_i) {
                result_type result = 0.0;
                for (size_t i = 0; i < m_k; i++) {
                        if (i != except_i) {
                                result += m_proposal[i];
                        }
                }
                return result;
        }
        template<class Engine>
        size_t draw_coordinate(Engine& eng, size_t except_i) {
                boost::random::uniform_int_distribution<> dist(0,m_k-2);
                size_t result = dist(eng);
                return result >= except_i ? result+1 : result;
        }
        template<class Engine>
        void draw_sample(Engine& eng, input_type sigma) {
                // initialize proposal distribution
                normal_distribution<input_type> dist(0.0, sigma);
                // copy the old state
                m_proposal = m_state;
                for (size_t i = 0; i < m_k; i++) {
                        // select the second coordinate at random
                        size_t j = draw_coordinate(eng, i);
                        // compute the range
                        result_type r = m_state[i] + m_state[j];
                        // draw a proposal
                        m_proposal[i] = (m_state[i] + r*dist(eng)) % r;
                        m_proposal[j] = 1.0 - sum_proposal(i);
                        // accept or reject
                        if (f(m_proposal) > 0.0 || )
                }
        }
        result_type f(const std::vector<result_type>& x) {
                result_type h = entropy(x);
                return m_beta(h)*histogram.pdf(h);
        }

        size_t m_k;
        input_type m_a1;
        input_type m_a2;
        std::vector<result_type> m_state;
        std::vector<result_type> m_proposal;
        boost::math::beta_distribution<input_type> m_beta;
        histogram_t<input_type, result_type> histogram;
};

} // namespace random
} // namespace boost


#endif /* __TFBAYES_ENTROPY_SAMPLER_HH__ */
