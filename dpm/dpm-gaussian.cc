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

DPM_GAUSSIAN::DPM_GAUSSIAN(
        gsl_matrix* Sigma,
        gsl_matrix* Sigma_0,
        gsl_vector* mu_0,
        const Data& data)
        : // raw data and cluster manager
          _data(data),
          _cluster_manager(data, new BivariateNormal(Sigma, Sigma_0, mu_0)),
          // strength parameter for the dirichlet process
          alpha(alpha)
{
}

DPM_GAUSSIAN::~DPM_GAUSSIAN() {
}

DPM_GAUSSIAN*
DPM_GAUSSIAN::clone() const {
        return new DPM_GAUSSIAN(*this);
}

bool
DPM_GAUSSIAN::valid_for_sampling(const element_t& element, const word_t& word) const
{
        return true;
}

void
DPM_GAUSSIAN::add_word(const word_t& word, cluster_tag_t tag)
{
        _cluster_manager[tag].add_word(word);
}

void
DPM_GAUSSIAN::remove_word(const word_t& word, cluster_tag_t tag)
{
        _cluster_manager[tag].remove_word(word);
}

size_t
DPM_GAUSSIAN::mixture_components() const
{
        return _cluster_manager.size() ;
}

void
DPM_GAUSSIAN::mixture_weights(const word_t& word, double weights[], cluster_tag_t tags[])
{
        size_t components = mixture_components();
        double sum        = -HUGE_VAL;
        double N = _data.size() - 1;

        cluster_tag_t i = 0;
        for (ClusterManager::const_iterator it = _cluster_manager.begin(); it != _cluster_manager.end(); it++) {
                Cluster& cluster = **it;
                tags[i] = cluster.tag();
                double num_elements = (double)cluster.size();
                weights[i] = log(num_elements/(alpha + N)) + cluster.distribution().log_pdf(word);
                // normalization constant
                sum = logadd(sum, weights[i]);
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        tags[components]    = _cluster_manager.get_free_cluster().tag();
        weights[components] = log(alpha/(alpha + N)) + _cluster_manager[tags[components]].distribution().log_pdf(word);
        sum = logadd(sum, weights[components]);

        ////////////////////////////////////////////////////////////////////////
        // normalize
        for (size_t i = 0; i < components+1; i++) {
                weights[i] = exp(weights[i] - sum);
        }
}

double
DPM_GAUSSIAN::likelihood() const {
        double result = 0;

        for (ClusterManager::const_iterator it = _cluster_manager.begin();
             it != _cluster_manager.end(); it++) {
                Cluster& cluster = **it;
                result += cluster.distribution().log_likelihood();
        }
        return result;
}

void
DPM_GAUSSIAN::update_posterior(size_t sampling_steps) {
}

const posterior_t&
DPM_GAUSSIAN::posterior() const {
        return _posterior;
}

const Data&
DPM_GAUSSIAN::data() const {
        return _data;
}

const ClusterManager&
DPM_GAUSSIAN::cluster_manager() const {
        return _cluster_manager;
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const DPM_GAUSSIAN& dpm)
{
        return o;
}
