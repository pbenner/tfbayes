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
#include <tfbayes/fastlog.h>

using namespace std;

DPM_Gaussian::DPM_Gaussian(
        double alpha,
        gsl_matrix* Sigma,
        gsl_matrix* Sigma_0,
        gsl_vector* mu_0,
        const DataGaussian& data)
        : _data(data),
          _cluster_assignments(_data.elements(), -1),
          _clustermanager(_cluster_assignments),
          // strength parameter for the dirichlet process
          alpha(alpha)
{
        _model_tag = _clustermanager.add_baseline_model(new BivariateNormal(Sigma, Sigma_0, mu_0, data));
        cluster_tag_t cluster_tag = _clustermanager.get_free_cluster(_model_tag).cluster_tag();
        for (DataGaussian::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                _clustermanager[cluster_tag].add_observations(range_t(*it, *it));
        }
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
DPM_Gaussian::add(const index_t& index, cluster_tag_t cluster_tag)
{
        _clustermanager[cluster_tag].add_observations(range_t(index,index));
}

void
DPM_Gaussian::remove(const index_t& index, cluster_tag_t cluster_tag)
{
        _clustermanager[cluster_tag].remove_observations(range_t(index,index));
}

size_t
DPM_Gaussian::mixture_components() const
{
        return _clustermanager.size();
}

void
DPM_Gaussian::mixture_weights(const index_t& index, double log_weights[], cluster_tag_t cluster_tags[])
{
        size_t components = mixture_components();
        double sum        = -HUGE_VAL;
        double N          = _data.elements() - 1;
        const range_t range(index, index);

        cluster_tag_t i = 0;
        for (ClusterManager::const_iterator it = _clustermanager.begin(); it != _clustermanager.end(); it++) {
                const Cluster& cluster = **it;
                cluster_tags[i] = cluster.cluster_tag();
                double num_elements = (double)cluster.size();
                // normalization constant
                sum = logadd(sum, log(num_elements/(alpha + N)) + cluster.model().log_pdf(range));
                log_weights[i] = sum;
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        cluster_tags[components] = _clustermanager.get_free_cluster(_model_tag).cluster_tag();
        sum = logadd(sum, log(alpha/(alpha + N)) + _clustermanager[cluster_tags[components]].model().log_pdf(range));
        log_weights[components] = sum;
}

gsl_matrix*
DPM_Gaussian::means() const {
        if (mixture_components() == 0) {
                return NULL;
        }

        gsl_matrix* means = gsl_matrix_alloc(mixture_components(), 2);

        size_t i = 0;
        for (ClusterManager::const_iterator it = _clustermanager.begin();
             it != _clustermanager.end(); it++) {
                const Cluster& cluster = **it;
                const BivariateNormal& bg = static_cast<const BivariateNormal&>(cluster.model());
                gsl_matrix_set(means, i, 0, gsl_vector_get(bg.mean(), 0));
                gsl_matrix_set(means, i, 1, gsl_vector_get(bg.mean(), 1));
                i++;
        }

        return means;
}

double
DPM_Gaussian::likelihood() const {
        double result = 0;

        for (ClusterManager::const_iterator it = _clustermanager.begin();
             it != _clustermanager.end(); it++) {
                Cluster& cluster = **it;
                result += cluster.model().log_likelihood();
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
DPM_Gaussian::clustermanager() const {
        return _clustermanager;
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const DPM_Gaussian& dpm)
{
        return o;
}
