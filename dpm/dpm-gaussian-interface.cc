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
#include <dpm-gaussian-interface.hh>
#include <dpm-tfbs.hh>
#include <sampler.hh>

using namespace std;

static DPM_TFBS* _gdpm;
static DataTFBS* _data;
static GibbsSampler* _sampler;

__BEGIN_C_REGION;

void _dpm_gaussian_init(double alpha, double lambda, int tfbs_length, int n, char *sequences[])
{
        __dpm_init__();

        vector<string> _sequences;
        for (size_t i = 0; i < (size_t)n; i++) {
                _sequences.push_back(sequences[i]);
        }

        _data    = new DataTFBS(_sequences);
        _gdpm    = new DPM_TFBS(alpha, lambda, (size_t)tfbs_length, *_data);
        _sampler = new GibbsSampler(*_gdpm, *_data);
}

unsigned int _dpm_gaussian_num_clusters() {
        return _gdpm->cluster_manager().size();
}

Bayes::Matrix* _dpm_gaussian_get_posterior() {
        Bayes::Matrix* result;
        const vector<vector<double> >& posterior = _gdpm->posterior();
        size_t n = posterior.size();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < posterior[i].size()) {
                        m = posterior[i].size();
                }
        }

        // allocate matrix
        result = Bayes::allocMatrix(n, m);
        // copy posterior
        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < m; j++) {
                        if (j < posterior[i].size()) {
                                result->mat[i][j] = posterior[i][j];
                        }
                        else {
                                result->mat[i][j] = 0;
                        }
                }
        }

        return result;
}

Bayes::Matrix* _dpm_gaussian_cluster_assignments() {
        Bayes::Matrix* result;
        size_t n = _gdpm->data().length();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < _gdpm->data().length(i)) {
                        m = _gdpm->data().length(i);
                }
        }

        // allocate matrix
        result = Bayes::allocMatrix(n, m);
        // copy posterior
        for (DataTFBS::const_iterator it = _gdpm->data().begin();
             it != _gdpm->data().end(); it++) {
                const index_t& index = *it;
                result->mat[index[0]][index[1]] = _gdpm->cluster_manager()[index];
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
