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
protected:
        // type definitions
        ////////////////////////////////////////////////////////////////////////
        typedef std::vector< input_type> ivector_t;
        typedef std::vector<result_type> rvector_t;
public:
        typedef entropy_distribution_t<input_type, result_type> edist_t;
protected:
        // member variables
        ////////////////////////////////////////////////////////////////////////
        edist_t m_entropy_distribution;
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
};

template <class input_type, class result_type, class counts_type>
result_type pdf(const entropy_multinomial_distribution_t<input_type, result_type>& dist,
                const std::vector<result_type>& theta,
                const counts_type& counts)
{
        assert(dist.k() == theta .size());
        assert(dist.k() == counts.size());

        // use static_cast to avoid constructing a temporary object
        const multinomial_distribution_t<result_type>& m =
                static_cast<const multinomial_distribution_t<result_type>&>(theta);

        return pdf(m, counts)*pdf(dist.entropy_distribution(), theta);
}


template <class input_type = double, class result_type = input_type>
class marginal_entropy_distribution_t : public entropy_multinomial_distribution_t<input_type, result_type>
{
        // type definitions
        ////////////////////////////////////////////////////////////////////////
        typedef entropy_multinomial_distribution_t<input_type, result_type> base_t;

        template <class T>
        class samples_cache_t : public std::vector<T>
        {
                typedef std::vector<T> base_t;
        public:
                samples_cache_t(size_t n)
                        : base_t(n, T())
                        { }
        };
public:
        typedef samples_cache_t<typename base_t::rvector_t> cache_t;
protected:
        // member variables
        ////////////////////////////////////////////////////////////////////////
        samples_cache_t<typename base_t::rvector_t> m_cache;
public:
        template <class Engine>
        marginal_entropy_distribution_t(size_t k, input_type a1, input_type a2, size_t n, Engine& eng)
                : base_t  (k, a1, a2)
                , m_cache (n) {

                for (typename cache_t::iterator it = m_cache.begin();
                     it != m_cache.end(); it++) {
                        *it = base_t::m_entropy_distribution(eng);
                }
        }
        const cache_t& cache() const {
                return m_cache;
        }
};

template <class input_type, class result_type, class counts_type>
result_type pdf(const marginal_entropy_distribution_t<input_type, result_type>& dist,
                const counts_type& counts)
{
        const entropy_multinomial_distribution_t<input_type, result_type>& base
                = static_cast<const entropy_multinomial_distribution_t<input_type, result_type>&>(dist);

        result_type result = 0.0;

        for (typename marginal_entropy_distribution_t<input_type, result_type>::cache_t::const_iterator it = dist.cache().begin();
             it != dist.cache().end(); it++) {
                result += pdf(base, *it, counts);
        }
        return result/dist.cache().size();
}

#endif /* __TFBAYES_ENTROPY_ENTROPY_MULTINOMIAL_DISTRIBUTION_HH__ */
