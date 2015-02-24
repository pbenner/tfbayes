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

        multinomial_distribution_t<result_type> m(theta);

        return pdf(m, counts)*pdf(dist.entropy_distribution(), theta);
}

template <class input_type, class result_type, class counts_type, class Engine>
result_type marginalize(entropy_multinomial_distribution_t<input_type, result_type>& dist,
                        const std::vector<counts_type>& counts,
                        Engine& eng, const size_t samples = 10000)
{
        std::vector<result_type> theta = dist.entropy_distribution()(eng);
        std::vector<result_type> proposal;
        // uniform distribution
        boost::random::uniform_01<input_type> runif;
        // marginal probability
        result_type result = 0.0;

        for (size_t i = 0; i < samples; i++) {
                proposal = dist.entropy_distribution()(eng);
                // accept or reject
                if (result_type(runif(eng)) <= std::min(
                            result_type(1.0),
                            pdf(dist, proposal, counts)/pdf(dist, theta, counts))) {
                        theta = proposal;
                }
        }
        for (size_t i = 0; i < samples; i++) {
                proposal = dist.entropy_distribution()(eng);
                // accept or reject
                if (result_type(runif(eng)) <= std::min(
                            result_type(1.0),
                            pdf(dist, proposal, counts)/pdf(dist, theta, counts))) {
                        theta = proposal;
                }
                // compute marginal
                result += pdf(dist, theta, counts);
        }
        return result/samples;
}

#endif /* __TFBAYES_ENTROPY_ENTROPY_MULTINOMIAL_DISTRIBUTION_HH__ */
