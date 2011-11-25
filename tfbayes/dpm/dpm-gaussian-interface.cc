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

#include <gsl/gsl_matrix.h>

namespace Bayes {
        extern "C" {
#include <tfbayes/linalg.h>
        }
}

#include <init.hh>
#include <data-gaussian.hh>
#include <dpm-gaussian-interface.hh>
#include <dpm-gaussian.hh>
#include <sampler.hh>

using namespace std;

static DPM_Gaussian* _gdpm;
static DataGaussian* _data;
static GibbsSampler* _sampler;

__BEGIN_C_REGION;

void _dpm_gaussian_init(
        int samples,
        double alpha,
        Bayes::Matrix* _Sigma,
        Bayes::Matrix* _Sigma_0,
        Bayes::Vector* _mu_0,
        Bayes::Vector* _pi)
{
        __dpm_init__();

        gsl_matrix *Sigma    = toGslMatrix(_Sigma);
        gsl_matrix *Sigma_0  = toGslMatrix(_Sigma_0);
        gsl_vector *mu_0     = toGslVector(_mu_0);
        const size_t cluster = _pi->size;

        _data    = new DataGaussian(cluster, (size_t)samples, Sigma, _pi->vec);
        _gdpm    = new DPM_Gaussian(alpha, Sigma, Sigma_0, mu_0, *_data);
        _sampler = new GibbsSampler(*_gdpm, *_data);

        gsl_matrix_free(Sigma);
        gsl_matrix_free(Sigma_0);
        gsl_vector_free(mu_0);
}

unsigned int _dpm_gaussian_num_clusters() {
        return _gdpm->clustermanager().size();
}

Bayes::Vector* _dpm_gaussian_cluster_tags() {
        size_t cluster = _gdpm->clustermanager().size();
        Bayes::Vector* cluster_tags = Bayes::allocVector(cluster);
        const ClusterManager& clustermanager = _gdpm->clustermanager();

        size_t i = 0;
        for (ClusterManager::const_iterator it = clustermanager.begin();
             it != clustermanager.end(); it++) {
                cluster_tags->vec[i++] = (*it)->cluster_tag();
        }
        return cluster_tags;
}

Bayes::Matrix* _dpm_gaussian_cluster_elements(int tag) {
        const ClusterManager& clustermanager = _gdpm->clustermanager();
        const Cluster& cluster = clustermanager[(cluster_tag_t)tag];
        Bayes::Matrix* elements = Bayes::allocMatrix(cluster.size(), 2);

        size_t i = 0;
        for (DataGaussian::const_iterator it = _data->begin();
             it != _data->end(); it++) {
                if (clustermanager[**it] == (cluster_tag_t)tag) {
                        elements->mat[i][0] = (*_data)[**it][0];
                        elements->mat[i][1] = (*_data)[**it][1];
                        i++;
                }
        }

        return elements;
}

Bayes::Matrix* _dpm_gaussian_get_posterior() {
        Bayes::Matrix* result;
        const vector<vector<double> >& probabilities = _gdpm->posterior().probabilities;
        size_t n = probabilities.size();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < probabilities[i].size()) {
                        m = probabilities[i].size();
                }
        }

        // allocate matrix
        result = Bayes::allocMatrix(n, m);
        // copy posterior
        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < m; j++) {
                        if (j < probabilities[i].size()) {
                                result->mat[i][j] = probabilities[i][j];
                        }
                        else {
                                result->mat[i][j] = 0;
                        }
                }
        }

        return result;
}

Bayes::Vector* _dpm_gaussian_cluster_assignments() {
        Bayes::Vector* result;
        size_t n = _data->elements();

        // allocate matrix
        result = Bayes::allocVector(n);
        // copy posterior
        for (DataGaussian::const_iterator it = _data->begin();
             it != _data->end(); it++) {
                const index_t& index = **it;
                result->vec[index[0]] = _gdpm->clustermanager()[index];
        }
        return result;
}

Bayes::Vector* _dpm_gaussian_hist_likelihood() {
        const vector<double>& likelihood = _sampler->sampling_history().likelihood[0];
        size_t length = likelihood.size();
        Bayes::Vector* result = Bayes::allocVector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = likelihood[i];
        }

        return result;
}

Bayes::Vector* _dpm_gaussian_hist_switches() {
        const vector<double>& switches = _sampler->sampling_history().switches[0];
        size_t length = switches.size();
        Bayes::Vector* result = Bayes::allocVector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = switches[i];
        }

        return result;
}

Bayes::Matrix* _dpm_gaussian_means() {
        gsl_matrix* tmp = _gdpm->means();
        Bayes::Matrix* means;
        if (tmp == NULL) {
                means = Bayes::allocMatrix(0,0);
        }
        else {  
                means = Bayes::fromGslMatrix(tmp);
                gsl_matrix_free(tmp);
        }
        return means;
}

Bayes::Matrix* _dpm_gaussian_data() {
        Bayes::Matrix* data = Bayes::allocMatrix(_data->elements(), 2);

        size_t i = 0;
        for (DataGaussian::const_iterator it = _data->begin();
             it != _data->end(); it++) {
                data->mat[i][0] = (*_data)[**it][0];
                data->mat[i][1] = (*_data)[**it][1];
                i++;
        }

        return data;
}

Bayes::Matrix* _dpm_gaussian_original_means() {
        gsl_matrix* means = _data->original_means();
        return Bayes::fromGslMatrix(means);
}

Bayes::Vector* _dpm_gaussian_original_cluster_assignments() {
        gsl_vector* original_cluster_assignments = _data->original_cluster_assignments();
        return Bayes::fromGslVector(original_cluster_assignments);
}

void _dpm_gaussian_print() {
        cout << *_gdpm << endl;
}

void _dpm_gaussian_sample(unsigned int n, unsigned int burnin) {
        _sampler->sample((size_t)n, (size_t)burnin);
}

void _dpm_gaussian_free() {
        delete(_data);
        delete(_gdpm);
        delete(_sampler);
}

__END_C_REGION;
