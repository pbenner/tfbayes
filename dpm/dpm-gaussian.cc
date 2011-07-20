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

#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <dpm-gaussian.hh>
#include <tfbayes/logarithmetic.h>

using namespace std;

DPM_Gaussian::DPM_Gaussian(
        gsl_matrix* Sigma,
        gsl_matrix* Sigma_0,
        gsl_vector* mu_0,
        const Data_Gaussian& data)
        : _data(data),
          _cluster_assignments(_data.length(), -1),
          _cluster_manager(new BivariateNormal(Sigma, Sigma_0, mu_0, data), _cluster_assignments),
          // strength parameter for the dirichlet process
          alpha(alpha)
{
}

DPM_Gaussian::~DPM_Gaussian() {
}

DPM_Gaussian*
DPM_Gaussian::clone() const {
        return new DPM_Gaussian(*this);
}

bool
DPM_Gaussian::valid_for_sampling(const index_t& index) const
{
        return true;
}

void
DPM_Gaussian::add(const index_t& index, cluster_tag_t tag)
{
        _cluster_manager[tag].add_observations(range_t(index,index));
}

void
DPM_Gaussian::remove(const index_t& index, cluster_tag_t tag)
{
        _cluster_manager[tag].remove_observations(range_t(index,index));
}

size_t
DPM_Gaussian::mixture_components() const
{
        return _cluster_manager.size();
}

void
DPM_Gaussian::mixture_weights(const index_t& index, double weights[], cluster_tag_t tags[])
{
        size_t components = mixture_components();
        double sum        = -HUGE_VAL;
        double N          = _data.size() - 1;
        const range_t range(index, index);

        cluster_tag_t i = 0;
        for (ClusterManager::const_iterator it = _cluster_manager.begin(); it != _cluster_manager.end(); it++) {
                Cluster& cluster = **it;
                tags[i] = cluster.tag();
                double num_elements = (double)cluster.size();
                weights[i] = log(num_elements/(alpha + N)) + cluster.distribution().log_pdf(range);
                // normalization constant
                sum = logadd(sum, weights[i]);
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        tags[components]    = _cluster_manager.get_free_cluster().tag();
        weights[components] = log(alpha/(alpha + N)) + _cluster_manager[tags[components]].distribution().log_pdf(range);
        sum = logadd(sum, weights[components]);

        ////////////////////////////////////////////////////////////////////////
        // normalize
        for (size_t i = 0; i < components+1; i++) {
                weights[i] = exp(weights[i] - sum);
        }
}

double
DPM_Gaussian::likelihood() const {
        double result = 0;

        for (ClusterManager::const_iterator it = _cluster_manager.begin();
             it != _cluster_manager.end(); it++) {
                Cluster& cluster = **it;
                result += cluster.distribution().log_likelihood();
        }
        return result;
}

void
DPM_Gaussian::update_posterior(size_t sampling_steps) {
}

const posterior_t&
DPM_Gaussian::posterior() const {
        return _posterior;
}

const ClusterManager&
DPM_Gaussian::cluster_manager() const {
        return _cluster_manager;
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const DPM_Gaussian& dpm)
{
        return o;
}
