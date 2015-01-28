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

#include <tfbayes/utility/histogram.hh>

namespace boost { namespace random {

template <class input_type = double, class result_type = input_type>
class entropy_distribution
{
public:
        entropy_distribution(size_t k, double a1, double a2) :
                m_k(k), m_a1(a1), m_a2(a2),
                m_state (k, 0.0),
                m_beta  (a1, a2) {
        }

        template<class Engine>
        const std::vector<result_type>& operator()(Engine& eng) {

                return m_state;
        }

private:
        size_t m_k;
        double m_a1;
        double m_a2;
        std::vector<result_type> m_state;
        boost::math::beta_distribution<input_type> m_beta;
};

} // namespace random
} // namespace boost


#endif /* __TFBAYES_ENTROPY_SAMPLER_HH__ */
