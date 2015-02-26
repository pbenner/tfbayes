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

#ifndef __TFBAYES_ENTROPY_ENTROPY_DISTRIBUTION_HH__
#define __TFBAYES_ENTROPY_ENTROPY_DISTRIBUTION_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <fstream>
#include <stdexcept>
#include <string>

#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/dirichlet.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/dirichlet_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_01.hpp>

#include <tfbayes/entropy/entropy.hh>
#include <tfbayes/entropy/entropy-approximation.hh>
#include <tfbayes/utility/histogram.hh>
#include <tfbayes/utility/boost-random-shuffle.hh>


template <class input_type, class result_type>
class entropy_distribution_t;

template <class input_type, class result_type>
result_type pdf(const entropy_distribution_t<input_type, result_type>& dist, const std::vector<result_type>& theta);

template <class input_type = double, class result_type = input_type>
class entropy_distribution_t
{
        // type definitions
        ////////////////////////////////////////////////////////////////////////
        typedef std::vector< input_type> ivector_t;
        typedef std::vector<result_type> rvector_t;
        typedef boost::math  ::beta_distribution<input_type> dbeta_t;
        typedef boost::math  ::dirichlet_distribution<input_type, result_type> ddirichlet_t;
        typedef boost::random::dirichlet_distribution<input_type, result_type> rdirichlet_t;
        // member variables
        ////////////////////////////////////////////////////////////////////////
        size_t m_k;
        // point in the probability simplex
        std::vector<result_type> m_theta;
        // transformed coordinates in the unit interval
        std::vector<result_type> m_phi;
        // same for drawing a proposal
        std::vector<result_type> m_proposal_theta;
        std::vector<result_type> m_proposal_phi;
        // prior density function
        boost::math::beta_distribution<input_type> m_dbeta;
        // weight function
        histogram_t<input_type, result_type> m_histogram;
        // is this the first sample?
        bool m_first_sample;
        // record some statistics
        double m_samples;
        double m_accepted;
        // sampling distributions
        boost::random::uniform_01<input_type> m_runif;
        ////////////////////////////////////////////////////////////////////////
public:
        entropy_distribution_t() :
                m_k              (0),
                m_first_sample   (true),
                m_samples        (0.0),
                m_accepted       (0.0)
                { }
        entropy_distribution_t(size_t k, input_type a1, input_type a2) :
                m_k              (k),
                m_theta          (k, 0.0),
                m_phi            (k, 1.0),
                m_proposal_theta (k, 0.0),
                m_proposal_phi   (k, 0.0),
                m_dbeta          (a1, a2),
                m_histogram      (entropy_approximation<input_type, result_type>(k)),
                m_first_sample   (true),
                m_samples        (0.0),
                m_accepted       (0.0)
                { }
        entropy_distribution_t(size_t k, input_type a1, input_type a2,
                               histogram_t<input_type, result_type> histogram) :
                m_k              (k),
                m_theta          (k, 0.0),
                m_phi            (k, 1.0),
                m_proposal_theta (k, 0.0),
                m_proposal_phi   (k, 0.0),
                m_dbeta          (a1, a2),
                m_histogram      (histogram),
                m_first_sample   (true),
                m_samples        (0.0),
                m_accepted       (0.0)
                { }

        template<class Engine>
        const rvector_t& operator()(Engine& eng) {
                if (m_first_sample) {
                        rdirichlet_t rdirichlet(ivector_t(m_k, 1.0));
                        m_theta = rdirichlet(eng);
                        transform_backward(m_phi, m_theta);
                        m_first_sample = false;
                }
                else {
                        draw_sample(eng);
                }

                return m_theta;
        }
        double acceptance_ratio() const {
                return m_accepted/m_samples;
        }
        const histogram_t<input_type, result_type>& histogram() const {
                return m_histogram;
        }
        const boost::math::beta_distribution<input_type>& beta() const {
                return m_dbeta;
        }
        const size_t& k() const {
                return m_k;
        }
protected:
        rvector_t& transform_forward(
                rvector_t& theta,
                rvector_t& phi) {
                result_type sum = 0.0;
                for (size_t i = 0; i < m_k; i++) {
                        theta[i] = phi[i]*(result_type(1.0) - sum);
                        sum     += theta[i];
                }
                return theta;
        }
        rvector_t& transform_backward(
                rvector_t& phi,
                rvector_t& theta) {
                result_type sum = 0.0;
                for (size_t i = 0; i < m_k-1; i++) {
                        phi[i] = theta[i]/(result_type(1.0) - sum);
                        sum   += theta[i];
                }
                return theta;
        }
        template<class Engine>
        size_t draw_coordinate(Engine& eng) {
                boost::random::uniform_int_distribution<> dist(0,m_k-2);
                return dist(eng);
        }
        template<class Engine>
        void draw_sample(Engine& eng) {
                // select the second coordinate at random
                for (size_t j = 0; j < m_k-1; j++) {
                        draw_sample(eng, j);
                }
                return draw_sample(eng, draw_coordinate(eng));
        }
        template<class Engine>
        void draw_sample(Engine& eng, size_t j) {
                // copy the old state
                m_proposal_theta = m_theta;
                m_proposal_phi   = m_phi;
                // proposal distribution
                rdirichlet_t rbeta(ivector_t({1.0, m_k-j-1.0}));
                // draw a proposal
                m_proposal_phi[j] = rbeta(eng)[0];
                transform_forward(m_proposal_theta, m_proposal_phi);
                // accept or reject
                if (result_type(m_runif(eng)) <= std::min(
                            result_type(1.0),
                            pdf(*this, m_proposal_theta)/pdf(*this, m_theta))) {
                         m_theta = m_proposal_theta;
                         m_phi   = m_proposal_phi;
                        // update statistics
                        m_accepted += 1.0;
                }
                m_samples += 1.0;
                boost::random::random_shuffle(m_theta.begin(), m_theta.end(), eng);
                transform_backward(m_phi, m_theta);
        }
        // template<class Engine>
        // void draw_sample(Engine& eng, input_type sigma) {
        //         // initialize proposal distribution
        //         boost::random::normal_distribution<input_type> rnorm(0.0, sigma);
        //         boost::random::uniform_01<input_type> runif;
        //         for (size_t i = 0; i < m_k; i++) {
        //                 // copy the old state
        //                m_proposal = m_state;
        //                 // select the second coordinate at random
        //                 size_t j = draw_coordinate(eng, i);
        //                 // compute the range
        //                 result_type r = m_state[i] + m_state[j];
        //                 // draw a proposal
        //                 m_proposal[i] = (m_state[i] + r*rnorm(eng)) % r;
        //                 m_proposal[j] = 1.0 - sum_proposal(j);
        //                 // accept or reject
        //                 if (static_cast<result_type>(runif(eng)) <= std::min(
        //                             static_cast<result_type>(1.0), pdf(*this, m_proposal)/pdf(*this, m_state))) {
        //                         m_state = m_proposal;
        //                         // update statistics
        //                         m_accepted += 1.0;
        //                 }
        //                 m_samples += 1.0;
        //         }
        // }
};

template <class input_type, class result_type>
result_type pdf(const entropy_distribution_t<input_type, result_type>& dist, const std::vector<result_type>& theta) {
        using boost::math::pdf;
        input_type x = input_type(entropy(theta))/std::log(dist.k());
        return pdf(dist.beta(), x)/pdf(dist.histogram(), x);
}

template <class input_type, class result_type>
result_type log_pdf(const entropy_distribution_t<input_type, result_type>& dist, const std::vector<result_type>& theta) {
        using boost::math::log_pdf;
        input_type x = input_type(entropy(theta))/std::log(dist.k());
        return std::log(pdf(dist.beta(), x)) - log_pdf(dist.histogram(), x);
}

#endif /* __TFBAYES_ENTROPY_ENTROPY_DISTRIBUTION_HH__ */
