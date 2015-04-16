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

#ifndef __TFBAYES_ENTROPY_APPROXIMATION_MIXTURE_HH__
#define __TFBAYES_ENTROPY_APPROXIMATION_MIXTURE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

template <class real_type = double, class probability_type = real_type>
class proposal_distribution_mixture_t
{
        typedef boost::random::dirichlet_distribution<real_type, probability_type> rdirichlet_t;
        typedef boost::math  ::dirichlet_distribution<real_type                  > ddirichlet_t;

        // cardinality
        size_t m_k;
        // number of components
        size_t m_size;
        // use exact summation
        bool m_extended_precision;
        // sampling distributions and densities
        std::vector<rdirichlet_t> m_rdirichlet;
        std::vector<ddirichlet_t> m_ddirichlet;

        real_type m_mean(real_type alpha) {
                alpha = std::max(alpha, real_type(1.0e-10));
                return boost::math::digamma(m_k*alpha + 1.0) - boost::math::digamma(alpha + 1.0);
        }
        real_type m_dmean(real_type alpha) {
                alpha = std::max(alpha, real_type(1.0e-10));
                return m_k*boost::math::trigamma(m_k*alpha + 1.0) - boost::math::trigamma(alpha + 1.0);
        }
        real_type m_sigma(real_type alpha) {
                alpha = std::max(alpha, real_type(1.0e-10));
                return std::sqrt((alpha+1.0)/(m_k*alpha+1.0)*boost::math::trigamma(alpha + 1.0)
                                 - boost::math::trigamma(m_k*alpha + 1.0));
        }
        real_type compute_alpha(real_type alpha, real_type target) {
                return newton<real_type>(boost::bind(&proposal_distribution_mixture_t::m_mean,  this, _1),
                                         boost::bind(&proposal_distribution_mixture_t::m_dmean, this, _1),
                                         alpha, target);
        }
public:
        // public type definitions
        struct parameters_t {
                double alpha_min;
                double alpha_max;
                double n_sigma;
        };
        // constructors
        proposal_distribution_mixture_t(size_t k, const parameters_t& parameters, bool extended_precision = false)
                : m_k(k), m_size(0.0), m_extended_precision(extended_precision) {

                real_type alpha_min = parameters.alpha_min;
                real_type alpha_max = parameters.alpha_max;

                for (real_type alpha = alpha_max; alpha >= alpha_min;) {
                        // verbose
                        std::cout << boost::format("Adding distribution at alpha = %f (with mean entropy %f)")
                                % alpha % m_mean(alpha) << std::endl;
                        // add distributions
                        m_rdirichlet.push_back(rdirichlet_t(k, alpha));
                        m_ddirichlet.push_back(ddirichlet_t(k, alpha));
                        // compute new alpha
                        alpha = m_mean(alpha) - parameters.n_sigma*m_sigma(alpha) > 0.0 ?
                                compute_alpha(alpha, m_mean(alpha) - parameters.n_sigma*m_sigma(alpha)) :
                                0.0;
                        // increase number of components
                        m_size += 1;
                }
                std::cout << "Done." << std::endl;
        }
        // operators
        template <class Engine>
        std::vector<probability_type> operator()(Engine& eng) {
                boost::random::uniform_int_distribution<> rint(0, m_size-1);
                return m_rdirichlet[rint(eng)](eng, m_extended_precision);
        }
        // access methods
        const rdirichlet_t& rdirichlet(size_t i) const {
                return m_rdirichlet[i];
        }
        const ddirichlet_t& ddirichlet(size_t i) const {
                return m_ddirichlet[i];
        }
        const size_t& size() const {
                return m_size;
        }
};

template <class real_type = double, class probability_type = real_type>
probability_type pdf(const proposal_distribution_mixture_t<real_type, probability_type>& d, const std::vector<probability_type>& x) {
        probability_type result = 0.0;
                
        for (size_t i = 0; i < d.size(); i++) {
                result += from_log_scale(boost::math::log_pdf(d.ddirichlet(i), x));
        }
        return result/probability_type(d.size());
}


#endif /* __TFBAYES_ENTROPY_APPROXIMATION_MIXTURE_HH__ */
