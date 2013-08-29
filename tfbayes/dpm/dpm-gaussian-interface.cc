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

#include <tfbayes/interface/common.hh>
#include <tfbayes/dpm/init.hh>
#include <tfbayes/dpm/data-gaussian.hh>
#include <tfbayes/dpm/dpm-gaussian.hh>
#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/utility/linalg.h>

using namespace std;

static dpm_gaussian_t* _gdpm;
static data_gaussian_t* _data;
static gibbs_sampler_t* _sampler;

__BEGIN_DECLS

// python interface
// -----------------------------------------------------------------------------

void _dpm_gaussian_init(
        int samples,
        double alpha,
        matrix_t* _Sigma,
        matrix_t* _Sigma_0,
        vector_t* _mu_0,
        vector_t* _pi)
{
        __dpm_init__();

        gsl_matrix *Sigma    = to_gsl_matrix(_Sigma);
        gsl_matrix *Sigma_0  = to_gsl_matrix(_Sigma_0);
        gsl_vector *mu_0     = to_gsl_vector(_mu_0);
        const size_t cluster = _pi->size;

        _data    = new data_gaussian_t(cluster, (size_t)samples, Sigma, _pi->vec);
        _gdpm    = new dpm_gaussian_t(alpha, Sigma, Sigma_0, mu_0, *_data);
        _sampler = new gibbs_sampler_t(*_gdpm, *_gdpm, *_data);

        gsl_matrix_free(Sigma);
        gsl_matrix_free(Sigma_0);
        gsl_vector_free(mu_0);
}

unsigned int _dpm_gaussian_num_clusters() {
        return _gdpm->state().size();
}

vector_t* _dpm_gaussian_cluster_tags() {
        size_t cluster = _gdpm->state().size();
        vector_t* cluster_tags = alloc_vector(cluster);
        const mixture_state_t& state = _gdpm->state();

        size_t i = 0;
        for (mixture_state_t::const_iterator it = state.begin();
             it != state.end(); it++) {
                cluster_tags->vec[i++] = (*it)->cluster_tag();
        }
        return cluster_tags;
}

matrix_t* _dpm_gaussian_cluster_elements(int tag) {
        const mixture_state_t& state = _gdpm->state();
        const cluster_t& cluster = state[(cluster_tag_t)tag];
        matrix_t* elements = alloc_matrix(cluster.size(), 2);

        size_t i = 0;
        for (data_gaussian_t::const_iterator it = _data->begin();
             it != _data->end(); it++) {
                if ((*_gdpm)[**it] == (cluster_tag_t)tag) {
                        elements->mat[i][0] = (*_data)[**it][0];
                        elements->mat[i][1] = (*_data)[**it][1];
                        i++;
                }
        }

        return elements;
}

matrix_t* _dpm_gaussian_get_posterior() {
        matrix_t* result;
        const vector<vector<double> >& probabilities = _gdpm->samples().probabilities;
        size_t n = probabilities.size();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < probabilities[i].size()) {
                        m = probabilities[i].size();
                }
        }

        // allocate matrix
        result = alloc_matrix(n, m);
        // copy samples
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

vector_t* _dpm_gaussian_cluster_assignments() {
        vector_t* result;
        size_t n = _data->elements();

        // allocate matrix
        result = alloc_vector(n);
        // copy samples
        for (data_gaussian_t::const_iterator it = _data->begin();
             it != _data->end(); it++) {
                const index_i& index = **it;
                result->vec[index[0]] = (*_gdpm)[index];
        }
        return result;
}

vector_t* _dpm_gaussian_hist_likelihood() {
        const vector<double>& likelihood = _sampler->sampling_history().likelihood[0];
        size_t length = likelihood.size();
        vector_t* result = alloc_vector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = likelihood[i];
        }

        return result;
}

vector_t* _dpm_gaussian_hist_switches() {
        const vector<double>& switches = _sampler->sampling_history().switches[0];
        size_t length = switches.size();
        vector_t* result = alloc_vector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = switches[i];
        }

        return result;
}

matrix_t* _dpm_gaussian_means() {
        gsl_matrix* tmp = _gdpm->means();
        matrix_t* means;
        if (tmp == NULL) {
                means = alloc_matrix(0,0);
        }
        else {  
                means = from_gsl_matrix(tmp);
                gsl_matrix_free(tmp);
        }
        return means;
}

matrix_t* _dpm_gaussian_data() {
        matrix_t* data = alloc_matrix(_data->elements(), 2);

        size_t i = 0;
        for (data_gaussian_t::const_iterator it = _data->begin();
             it != _data->end(); it++) {
                data->mat[i][0] = (*_data)[**it][0];
                data->mat[i][1] = (*_data)[**it][1];
                i++;
        }

        return data;
}

matrix_t* _dpm_gaussian_original_means() {
        gsl_matrix* means = _data->original_means();
        return from_gsl_matrix(means);
}

vector_t* _dpm_gaussian_original_cluster_assignments() {
        gsl_vector* original_cluster_assignments = _data->original_cluster_assignments();
        return from_gsl_vector(original_cluster_assignments);
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

__END_DECLS
