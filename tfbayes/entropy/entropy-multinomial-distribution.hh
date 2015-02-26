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

#ifndef __TFBAYES_ENTROPY_ENTROPY_MULTINOMIAL_DISTRIBUTION_HH__
#define __TFBAYES_ENTROPY_ENTROPY_MULTINOMIAL_DISTRIBUTION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <boost/random/uniform_01.hpp>

#include <tfbayes/entropy/entropy-distribution.hh>
#include <tfbayes/entropy/multinomial-distribution.hh>
#include <tfbayes/utility/linalg.hh>

template <class input_type = double, class result_type = input_type>
class entropy_multinomial_distribution_t
{
        // type definitions
        ////////////////////////////////////////////////////////////////////////
        typedef std::vector< input_type> ivector_t;
        typedef std::vector<result_type> rvector_t;
public:
        typedef entropy_distribution_t<input_type, result_type> edist_t;
protected:
        // member variables
        ////////////////////////////////////////////////////////////////////////
        entropy_distribution_t<input_type, result_type> m_entropy_distribution;
public:
        entropy_multinomial_distribution_t(size_t k, input_type a1, input_type a2)
                : m_entropy_distribution (k, a1, a2)
                { }

        const edist_t& entropy_distribution() const {
                return m_entropy_distribution;
        }
        edist_t& entropy_distribution() {
                return m_entropy_distribution;
        }
        const size_t& k() const {
                return m_entropy_distribution.k();
        }
protected:
};

template <class input_type, class result_type, class counts_type>
result_type pdf(const entropy_multinomial_distribution_t<input_type, result_type>& dist,
                const std::vector<result_type>& theta,
                const std::vector<counts_type>& counts)
{
        assert(dist.k() == theta .size());
        assert(dist.k() == counts.size());

        // use static_cast to avoid constructing a temporary object
        const multinomial_distribution_t<result_type>& m =
                static_cast<const multinomial_distribution_t<result_type>&>(theta);

        return pdf(m, counts)*pdf(dist.entropy_distribution(), theta);
}

template <class input_type, class result_type, class counts_type, class Engine>
result_type marginalize(entropy_multinomial_distribution_t<input_type, result_type>& dist,
                        const std::vector<counts_type>& counts,
                        size_t samples, Engine& eng)
{
        result_type result = 0.0;

        for (size_t i = 0; i < samples; i++) {
                result += pdf(dist, dist.entropy_distribution()(eng), counts);
        }
        return result/samples;
}

template <class input_type, class result_type, class counts_type, class cache_type>
result_type marginalize(const entropy_multinomial_distribution_t<input_type, result_type>& dist,
                        const std::vector<counts_type>& counts,
                        const cache_type& samples_cache)
{
        result_type result = 0.0;

        for (typename cache_type::const_iterator it = samples_cache.begin();
             it != samples_cache.end(); it++) {
                result += pdf(dist, *it, counts);
        }
        return result/samples_cache.size();
}

template <class input_type, class result_type, class cache_type, class Engine>
void marginalize_fill_cache(entropy_multinomial_distribution_t<input_type, result_type>& dist,
                            cache_type& samples_cache, Engine& eng)
{
        for (typename cache_type::iterator it = samples_cache.begin();
             it != samples_cache.end(); it++) {
                *it = dist.entropy_distribution()(eng);
        }
}

template <class T>
class samples_cache_t : public std::vector<T>
{
        typedef std::vector<T> base_t;
public:
        samples_cache_t(size_t n)
                : base_t(n, T())
                { }
};

#endif /* __TFBAYES_ENTROPY_ENTROPY_MULTINOMIAL_DISTRIBUTION_HH__ */
