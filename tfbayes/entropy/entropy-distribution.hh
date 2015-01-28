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

#include <tfbayes/utility/histogram.hh>

namespace boost { namespace random {

template <class input_type = double, class result_type = input_type>
class entropy_distribution
{
public:
        entropy_distribution(size_t k, double a1, double a2) :
                m_size(alpha.size()) {
                BOOST_FOREACH(const input_type& a, alpha) {
                        m_distributions.push_back(
                                gamma_distribution<input_type>(a));
                }
        }

        template<class Engine>
        std::vector<result_type> operator()(Engine& eng) {
        }

private:
        size_t m_size;
        std::vector<double> m_state;
};

} // namespace random
} // namespace boost


#endif /* __TFBAYES_ENTROPY_SAMPLER_HH__ */
