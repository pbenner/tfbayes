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

#ifndef __TFBAYES_ENTROPY_MULTINOMIAL_DISTRIBUTION_HH__
#define __TFBAYES_ENTROPY_MULTINOMIAL_DISTRIBUTION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cmath>


template <class result_type = double>
class multinomial_distribution_t : public std::vector<result_type>
{
        // derive multinomial_distribution_t from std::vector to allow
        // static_casts
public:
        typedef std::vector<result_type> base_t;

        multinomial_distribution_t(const base_t& theta)
                : base_t(theta)
                { }
        template <class T>
        multinomial_distribution_t(const T& theta)
                : base_t(theta.begin(), theta.end())
                { }

        multinomial_distribution_t(const multinomial_distribution_t& m)
                : base_t(m.begin(), m.end())
                { }

        size_t k() const {
                return base_t::size();
        }
        const base_t& theta() const {
                return *this;
        }
};

template <class result_type, class T>
result_type pdf(const multinomial_distribution_t<result_type>& dist, const T& counts)
{
        assert(dist.k() == counts.size());

        result_type result = 1.0;

        for (size_t i = 0; i < dist.k(); i++) {
                result *= std::pow(dist.theta()[i], counts[i]);
        }
        return result;
}

#endif /* __TFBAYES_ENTROPY_MULTINOMIAL_DISTRIBUTION_HH__ */
