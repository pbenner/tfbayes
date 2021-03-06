/* Copyright (C) 2011 Philipp Benner
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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cmath>

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <tfbayes/dpm/dpm-tfbs-prior.hh>

using namespace std;

// Pitman-Yor prior
////////////////////////////////////////////////////////////////////////////////

pitman_yor_prior::pitman_yor_prior(
        double alpha, double discount)
        : dpm_tfbs_prior_t(),
          alpha(alpha),
          discount(discount)
{}

pitman_yor_prior*
pitman_yor_prior::clone() const {
        return new pitman_yor_prior(*this);
}

double
pitman_yor_prior::log_predictive(const cluster_t& cluster, const dpm_tfbs_state_t& state, size_t n) const
{
        const double N = state.num_tfbs;
        const double K = state.size() - state.bg_cluster_tags.size();
        double result = 0.0;

        for (size_t i = 0; i < n; i++) {
                if (cluster.size() + i == 0) {
                        result += log(alpha + discount*K) - log(N + i + alpha);
                }
                else {
                        result += log(cluster.size() + i - discount) - log(N + i + alpha);
                }
        }
        return result;
}

double
pitman_yor_prior::joint(const dpm_tfbs_state_t& state) const
{
        const double N = state.num_tfbs;
        const double K = state.size() - state.bg_cluster_tags.size();

        double sum = K*log(alpha) + boost::math::lgamma<double>(alpha) - boost::math::lgamma<double>(N + alpha);

        for (cl_iterator it = state.begin(); it != state.end(); it++) {
                const cluster_t& cluster = **it;
                if (!state.is_background(cluster)) {
                        sum += boost::math::lgamma<double>(cluster.size());
                }
        }

        return sum;
}

// Uniform prior
////////////////////////////////////////////////////////////////////////////////

uniform_prior::uniform_prior(double alpha)
        : dpm_tfbs_prior_t(),
          alpha(alpha)
{}

uniform_prior*
uniform_prior::clone() const {
        return new uniform_prior(*this);
}

double
uniform_prior::log_predictive(const cluster_t& cluster, const dpm_tfbs_state_t& state, size_t n) const
{
        assert(n == 1);
        const double K = state.size() - state.bg_cluster_tags.size();

        if (cluster.size() == 0) {
                return log(alpha) - log(alpha + K);
        }
        else {
                return -log(alpha + K);
        }
}

double
uniform_prior::joint(const dpm_tfbs_state_t& state) const
{
        return 0;
}

// Poppe prior
////////////////////////////////////////////////////////////////////////////////

poppe_prior::poppe_prior()
        : dpm_tfbs_prior_t()
{}

poppe_prior*
poppe_prior::clone() const {
        return new poppe_prior(*this);
}

double
poppe_prior::log_predictive(const cluster_t& cluster, const dpm_tfbs_state_t& state, size_t n) const
{
        assert(n == 1);
        const double K = state.size() - state.bg_cluster_tags.size();
        const double N = state.num_tfbs;

        if (K == 0 && cluster.size() == 0) {
                return 0;
        }
        if (cluster.size() == 0) {
                if (K == 1.0) {
                        return -log(N+1.0);
                }
                else {
                        return log(K*(K-1.0)/(N*(N+1.0)));
                }
        }
        else {
                if (K == 1.0) {
                        return log(N/(N+1.0));
                }
                else {
                        return log((cluster.size()+1.0)/(N+1.0) * (N-K+1.0)/N);
                }
        }
}

double
poppe_prior::joint(const dpm_tfbs_state_t& state) const
{
        return 0;
}
