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

#include <gsl/gsl_matrix.h>

namespace Bayes {
        extern "C" {
#include <tfbayes/linalg.h>
        }
}

#include <init.hh>
#include <interface.hh>
#include <sampler.hh>

using namespace std;

static DPM* _gdpm;
static Data* _data;
static GibbsSampler* _sampler;

__BEGIN_C_REGION;

Bayes::Vector * _allocVector(int size)              { return Bayes::allocVector(size); }
void            _freeVector(Bayes::Vector *v)       { Bayes::freeVector(v); }
Bayes::Matrix * _allocMatrix(int rows, int columns) { return Bayes::allocMatrix(rows, columns); }
void            _freeMatrix(Bayes::Matrix *m)       { Bayes::freeMatrix(m); }
void            _free(void *ptr)                    { free(ptr); }

void _dpm_init(int n, char *sequences[])
{
        __dpm_init__();

        _data = new Data(n, sequences);
        _gdpm = new DPM(*_data);
        _sampler = new GibbsSampler(*_gdpm, *_data);
}

unsigned int _dpm_num_clusters() {
//        return _gdpm->num_clusters();
        return 0;
}

Bayes::Matrix* _dpm_get_posterior() {
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

Bayes::Matrix* _dpm_cluster_assignments() {
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
        for (Data::const_iterator it = _gdpm->data().begin();
             it != _gdpm->data().end(); it++) {
                const element_t& element = *it;
                result->mat[element.sequence][element.position] = _gdpm->cluster_manager()[element];
        }
        return result;
}

Bayes::Vector* _dpm_hist_likelihood() {
        vector<double>& likelihood = _sampler->get_hist_likelihood();
        size_t length = likelihood.size();
        Bayes::Vector* result = Bayes::allocVector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = likelihood[i];
        }

        return result;
}

Bayes::Vector* _dpm_hist_switches() {
        vector<double>& switches = _sampler->get_hist_switches();
        size_t length = switches.size();
        Bayes::Vector* result = Bayes::allocVector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = switches[i];
        }

        return result;
}

void _dpm_print() {
//        cout << *_gdpm << endl;
}

void _dpm_sample(unsigned int n, unsigned int burnin) {
        _sampler->sample((size_t)n, (size_t)burnin);
}

void _dpm_free() {
        delete(_gdpm);
}

__END_C_REGION;
