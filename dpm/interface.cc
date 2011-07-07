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

#include "init.hh"
#include "clusters.hh"
#include "data.hh"
#include "dpm.hh"
#include "interface.hh"

static DPM* _gdpm;

__BEGIN_C_REGION;

Bayes::Vector * _allocVector(int size)              { return Bayes::allocVector(size); }
void            _freeVector(Bayes::Vector *v)       { Bayes::freeVector(v); }
Bayes::Matrix * _allocMatrix(int rows, int columns) { return Bayes::allocMatrix(rows, columns); }
void            _freeMatrix(Bayes::Matrix *m)       { Bayes::freeMatrix(m); }
void            _free(void *ptr)                    { free(ptr); }

void _dpm_init(int n, char *sequences[])
{
        __dpm_init__();

        _gdpm = new DPM((size_t)n, sequences);
}

unsigned int _dpm_num_clusters() {
//        return _gdpm->num_clusters();
        return 0;
}

Bayes::Matrix* _dpm_get_posterior() {
        Bayes::Matrix* result;
        const vector<vector<double> >& posterior = _gdpm->get_posterior();
        size_t n = posterior.size();
        size_t length = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (length < posterior[i].size()) {
                        length = posterior[i].size();
                }
        }

        // allocate matrix
        result = Bayes::allocMatrix(n, m);
        // copy posterior
        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < length; j++) {
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

Bayes::Matrix* _dpm_cluster(unsigned int c) {
        // Clusters::iterator it = (*_gdpm).get_clusters().begin();
        // for (size_t i = 0; i < c; i++) {
        //         it++;
        // }
        // Cluster::cluster& cl = **it;
        // int n = cl.elements.size();
        // int m = 2;

        // Bayes::Matrix* result = Bayes::allocMatrix(n, m);
        // int i = 0;
        // for (Clusters::elements_t::iterator it = cl.elements.begin();
        //      it != cl.elements.end(); it++) {
        //         result->mat[i][0] = (*it)->x[0];
        //         result->mat[i][1] = (*it)->x[1];
        //         i++;
        // }

        // return result;
        return NULL;
}

Bayes::Vector* _dpm_hist_likelihood() {
        vector<double>& likelihood = _gdpm->get_hist_likelihood();
        size_t length = likelihood.size();
        Bayes::Vector* result = Bayes::allocVector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = likelihood[i];
        }

        return result;
}

Bayes::Vector* _dpm_hist_switches() {
        vector<double>& switches = _gdpm->get_hist_switches();
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
        _gdpm->gibbs_sample((size_t)n, (size_t)burnin);
}

void _dpm_free() {
        delete(_gdpm);
}

__END_C_REGION;
