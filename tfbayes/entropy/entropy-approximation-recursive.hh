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

#ifndef __TFBAYES_ENTROPY_APPROXIMATION_RECURSIVE_HH__
#define __TFBAYES_ENTROPY_APPROXIMATION_RECURSIVE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/entropy/entropy-approximation.hh>
#include <tfbayes/entropy/entropy-distribution.hh>


template <class real_type = double, class probability_type = real_type>
class proposal_distribution_recursive_t
{
        // entropy distribution with cardinality m_k and a histogram
        // from the m_k-1 approximation
        entropy_distribution_t<real_type, probability_type> m_entropy_distribution;
public:
        proposal_distribution_recursive_t(size_t k)
                : m_entropy_distribution(k, 1, 1, entropy_approximation<real_type, probability_type>(k-1)) {
        }

        template <class Engine>
        std::vector<probability_type> operator()(Engine& eng) {
                return m_entropy_distribution(eng);
        }
        entropy_distribution_t<real_type, probability_type>& entropy_distribution() {
                return m_entropy_distribution;
        }
        const entropy_distribution_t<real_type, probability_type>& entropy_distribution() const {
                return m_entropy_distribution;
        }
};

template <class real_type = double, class probability_type = real_type>
probability_type pdf(const proposal_distribution_recursive_t<real_type, probability_type>& d, const std::vector<probability_type>& x) {
        return pdf(d.entropy_distribution(), x);
}

#endif /* __TFBAYES_ENTROPY_APPROXIMATION_RECURSIVE_HH__ */
