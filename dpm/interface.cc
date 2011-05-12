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

void _dpm_init(int n, int m, char *sequences[], int *clusters[])
{
        __dpm_init__();

        Data* data = new Data(n, m, sequences, clusters);
        _gdpm = new DPM(data);
}

unsigned int _dpm_num_clusters() {
        return _gdpm->num_clusters();
}

Bayes::Matrix* _dpm_get_posterior() {
        return Bayes::fromGslMatrix((*_gdpm).get_posterior());
}

Bayes::Matrix* _dpm_cluster(unsigned int c) {
        Clusters::iterator it = (*_gdpm).get_clusters().begin();
        for (unsigned int i = 0; i < c; i++) {
                it++;
        }
        Clusters::cluster& cl = **it;
        int n = cl.elements.size();
        int m = 2;

        Bayes::Matrix* result = Bayes::allocMatrix(n, m);
        int i = 0;
        for (Clusters::elements_t::iterator it = cl.elements.begin();
             it != cl.elements.end(); it++) {
                result->mat[i][0] = (*it)->x[0];
                result->mat[i][1] = (*it)->x[1];
                i++;
        }

        return result;
}

Bayes::Vector* _dpm_original_tags(unsigned int c) {
        Clusters::cluster& cl = (*_gdpm)[c];
        int n = cl.elements.size();

        Bayes::Vector* result = Bayes::allocVector(n);
        int i = 0;
        for (Clusters::elements_t::iterator it = cl.elements.begin();
             it != cl.elements.end(); it++) {
                result->vec[i] = (double)((*it)->original_cluster);
                i++;
        }

        return result;
}

Bayes::Vector* _dpm_hist_likelihood() {
        vector<double>& likelihood = _gdpm->get_hist_likelihood();
        unsigned int len = likelihood.size();
        Bayes::Vector* result = Bayes::allocVector(len);

        for (unsigned int i = 0; i < len; i++) {
                result->vec[i] = likelihood[i];
        }

        return result;
}

Bayes::Vector* _dpm_hist_switches() {
        vector<double>& switches = _gdpm->get_hist_switches();
        unsigned int len = switches.size();
        Bayes::Vector* result = Bayes::allocVector(len);

        for (unsigned int i = 0; i < len; i++) {
                result->vec[i] = switches[i];
        }

        return result;
}

void _dpm_print() {
        cout << *_gdpm << endl;
}

void _dpm_sample(unsigned int n) {
        _gdpm->gibbsSample(n);
}

void _dpm_free() {
        delete(_gdpm);
}

__END_C_REGION;
