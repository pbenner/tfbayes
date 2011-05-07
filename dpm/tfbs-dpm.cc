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

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "tfbs-dpm.hh"

using namespace std;

TfbsDPM::TfbsDPM(TfbsData* data)
        : DPM(data), lambda(0.99),
          pd_tfbs_alpha(gsl_matrix_alloc(10, 4)),
          pd_bg_alpha(gsl_matrix_alloc(1, 4))
{
        alpha = 0.01;

        for (int i = 0; i < 10; i++) {
                for (int j = 0; j < 4; j++) {
                        gsl_matrix_set(pd_tfbs_alpha, i, j, 1);
                }
        }
        for (int j = 0; j < 4; j++) {
                gsl_matrix_set(pd_bg_alpha, 0, j, 10);
        }

        // initialize distributions
        predictiveDist_tfbs          = new ProductDirichlet(lambda, pd_tfbs_alpha);
        posteriorPredictiveDist_tfbs = new ProductDirichlet();
        posteriorPredictiveDist_bg   = new ProductDirichlet();
}

TfbsDPM::~TfbsDPM() {
        delete(predictiveDist);
        delete(posteriorPredictiveDist);
}

void
TfbsDPM::count_statistic(const Cluster::cluster& cluster, gsl_matrix* alpha, gsl_matrix* counts) {
        int len = counts->size1;

        // reset counts
        for (int i = 0; i < len; i++) {
                gsl_matrix_set(counts, i, 0, gsl_matrix_get(alpha, i, 0));
                gsl_matrix_set(counts, i, 1, gsl_matrix_get(alpha, i, 1));
                gsl_matrix_set(counts, i, 2, gsl_matrix_get(alpha, i, 2));
                gsl_matrix_set(counts, i, 3, gsl_matrix_get(alpha, i, 3));
        }
        // compute count statistic
        for (Cluster::elements_t::const_iterator it  = cluster.elements.begin();
             it != cluster.elements.end(); it++) {
                char buf[len];
                ((TfbsData*)da)->get_nucleotide(**it, len, buf);
                for (int i = 0; i < len; i++) {
                        switch (buf[i]) {
                        case 'A':
                        case 'a':
                                gsl_matrix_set(counts, i, 0,
                                               gsl_matrix_get(counts, i, 0)+1);
                        break;
                        case 'C':
                        case 'c':
                                gsl_matrix_set(counts, i, 1,
                                               gsl_matrix_get(counts, i, 1)+1);
                        break;
                        case 'G':
                        case 'g':
                                gsl_matrix_set(counts, i, 2,
                                               gsl_matrix_get(counts, i, 2)+1);
                        break;
                        case 'T':
                        case 't':
                                gsl_matrix_set(counts, i, 3,
                                               gsl_matrix_get(counts, i, 3)+1);
                        break;
                        }
                }
        }
}

Distribution&
TfbsDPM::posteriorPredictive(const Cluster::cluster& cluster) {
        if (cluster.tag == 0) {
                // background model
                gsl_matrix* counts = gsl_matrix_alloc(1, 4);
                count_statistic(cluster, pd_bg_alpha, counts);
                ((ProductDirichlet *)posteriorPredictiveDist_bg)->update(1-lambda, counts);
                return *posteriorPredictiveDist_bg;
        }
        else {
                // motif model
                gsl_matrix* counts = gsl_matrix_alloc(10, 4);
                count_statistic(cluster, pd_tfbs_alpha, counts);
                ((ProductDirichlet *)posteriorPredictiveDist_tfbs)->update(lambda, counts);
                return *posteriorPredictiveDist_tfbs;
        }
}

Distribution&
TfbsDPM::predictive() {
        return *predictiveDist_tfbs;
}

double
TfbsDPM::likelihood() {
        return 0.0;
}

void
TfbsDPM::compute_statistics() {
}
