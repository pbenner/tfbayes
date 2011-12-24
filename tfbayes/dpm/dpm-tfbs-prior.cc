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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_sf_gamma.h>

#include <dpm-tfbs-prior.hh>

// Pitman-Yor prior
////////////////////////////////////////////////////////////////////////////////

pitman_yor_prior::pitman_yor_prior(
        const dpm_tfbs_state_t& state,
        double alpha, double discount)
        : state(state), alpha(alpha), discount(discount)
{}

pitman_yor_prior*
pitman_yor_prior::clone() const {
        return new pitman_yor_prior(*this);
}

double
pitman_yor_prior::predictive(const cluster_t& cluster) const
{
        const double K = state.size()-1;

        if (cluster.size() == 0) {
                return log(alpha + discount*(K)) - log(state.num_tfbs + alpha);
        }
        else {
                return log(cluster.size()-discount) - log(state.num_tfbs + alpha);
        }
}

double
pitman_yor_prior::likelihood() const
{
        const double N = state.num_tfbs;

        double sum = gsl_sf_lngamma(alpha) - gsl_sf_lngamma(N + alpha);

        for (cl_iterator it = state.begin(); it != state.end(); it++) {
                cluster_t& cluster = **it;
                sum += gsl_sf_lngamma(cluster.size());
        }

        return sum;
}

// Uniform prior
////////////////////////////////////////////////////////////////////////////////

uniform_prior::uniform_prior(
        const dpm_tfbs_state_t& state, double alpha)
        : state(state), alpha(alpha)
{}

uniform_prior*
uniform_prior::clone() const {
        return new uniform_prior(*this);
}

double
uniform_prior::predictive(const cluster_t& cluster) const
{
        const double K = state.size()-1;

        if (cluster.size() == 0) {
                return log(alpha) - log(alpha + K);
        }
        else {
                return -log(alpha + K);
        }
}

double
uniform_prior::likelihood() const
{
        return 0;
}

// Poppe prior
////////////////////////////////////////////////////////////////////////////////

poppe_prior::poppe_prior(
        const dpm_tfbs_state_t& state)
        : state(state)
{}

poppe_prior*
poppe_prior::clone() const {
        return new poppe_prior(*this);
}

double
poppe_prior::predictive(const cluster_t& cluster) const
{
        const double K = state.size()-1;
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
poppe_prior::likelihood() const
{
        return 0;
}
