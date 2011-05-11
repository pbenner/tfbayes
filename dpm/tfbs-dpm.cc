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

#include <tfbayes/logarithmetic.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "tfbs-dpm.hh"

using namespace std;

TfbsDPM::TfbsDPM(TfbsData* data)
        : DPM(data), lambda(0.4),
          pd_tfbs_alpha(gsl_matrix_alloc(10, 4)),
          pd_bg_alpha(gsl_matrix_alloc(1, 4))
{
        alpha = 1;

        for (int i = 0; i < 10; i++) {
                for (int j = 0; j < 4; j++) {
                        gsl_matrix_set(pd_tfbs_alpha, i, j, 1);
                }
        }
        for (int j = 0; j < 4; j++) {
                gsl_matrix_set(pd_bg_alpha, 0, j, 10);
        }

        // initialize distributions
        predictiveDist_tfbs          = new ProductDirichlet(pd_tfbs_alpha);
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
                ((ProductDirichlet *)posteriorPredictiveDist_bg)->update(counts);
                return *posteriorPredictiveDist_bg;
        }
        else {
                // motif model
                gsl_matrix* counts = gsl_matrix_alloc(10, 4);
                count_statistic(cluster, pd_tfbs_alpha, counts);
                ((ProductDirichlet *)posteriorPredictiveDist_tfbs)->update(counts);
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

bool
TfbsDPM::check_element(Data::element& element)
{
        // check if there is enough space for the binding site
        if (((TfbsData*)da)->num_successors(element) < 9) {
                return false;
        }
        if (cl.getClusterTag(element) > 0) {
                // this is alreade a binding site, continue
                return true;
        }
        if (cl.getClusterTag(element) == 0) {
                Data::iterator it = da->find(element);
                // check if all successing nucleotides belong to the background
                for (int i = 0; i < 9 && it != da->end(); i++) {
                        it++;
//                        printf("neighbours: %d\n", cl.getClusterTag(*it));
                        if (cl.getClusterTag(*it) != 0) {
                                return false;
                        }
                }
                // check if all previous nucleotides belong to the background
                it = da->find(element);
                for (int i = 0; i < 9 && it != da->begin(); i++) {
                        it--;
//                        printf("neighbours: %d\n", cl.getClusterTag(*it));
                        if (cl.getClusterTag(*it) != 0) {
                                return false;
                        }
                }
                // there is enough space to the left and right, continue
                return true;
        }
        // this element is within a binding site, abort
        return false;
}

void
TfbsDPM::release_block(Data::element& element) {
        Data::iterator it = da->find(element);
        for (int i = 0; i < 10; i++) {
                if (cl.getClusterTag(*it) != -1) {
                        cl.release(*it);
                }
                it++;
        }
}

void
TfbsDPM::assign_block(Data::element& element, Cluster::cluster_tag_t c) {
        if (c == 0) {
                Data::iterator it = da->find(element);
                for (int i = 0; i < 10; i++) {
                        cl.assign(*it, 0);
                        it++;
                }
        }
        else {
                cl.assign(element, c);
        }
}

bool
TfbsDPM::sample(Data::element& element) {
        if (check_element(element) == false) {
                return false;
        }

        Cluster::cluster_tag_t old_cluster_tag = cl.getClusterTag(element);
        release_block(element);
        Distribution& pred = predictive();
        Cluster::size_type num_clusters = cl.size();
        double weights[num_clusters+1];
        Cluster::cluster_tag_t tags[num_clusters+1];
        double sum = 0;

        char buf[10];
        ((TfbsData*)da)->get_nucleotide(element, 10, buf);

        ////////////////////////////////////////////////////////////////////////
        // mixture component 1: background model
        Cluster::iterator it = cl.begin();
        {
                Distribution& postPred = posteriorPredictive(**it);
                weights[0] = 1;
                for (int i = 0; i < 10; i++) {
                        weights[0] *= postPred.pdf(buf+i);
                }
                weights[0] *= (1-lambda);
                tags[0]     = (*it)->tag;
                sum = weights[0];
//                printf("_weights[%ld]: %.10f\n", 0, weights[0]);
        }

        ////////////////////////////////////////////////////////////////////////
        // mixture component 2: dirichlet process for tfbs models
        it++;
        for (Cluster::cluster_tag_t i = 1; it != cl.end(); i++) {
                Distribution& postPred = posteriorPredictive(**it);
                double num_elements    = (double)(*it)->elements.size();
                weights[i] = lambda*num_elements*postPred.pdf(buf);
                tags[i]    = (*it)->tag;
                // normalization constant
                sum += weights[i];
                it++;
        }
        // add the tag of a new class and compute their weight
//        printf("alpha: %f\n", alpha);
        weights[num_clusters] = alpha*pred.pdf(buf);
        tags[num_clusters]    = cl.next_free_cluster()->tag;
        sum += weights[num_clusters];

        // normalize
//        printf("sum: %.10f\n", sum);
        for (Cluster::cluster_tag_t i = 0; i < num_clusters+1; i++) {
//                printf("weights[%ld]: %.10f\n", i, weights[i]);
                weights[i] /=  sum;
//                printf("weights[%ld]: %.10f\n", i, weights[i]);
        }

        // draw a new cluster for the element
        gsl_ran_discrete_t* gdd  = gsl_ran_discrete_preproc(num_clusters+1, weights);
        Cluster::cluster_tag_t i = gsl_ran_discrete(_r, gdd);
        gsl_ran_discrete_free(gdd);

        assign_block(element, tags[i]);
        printf("\n\n");

        return old_cluster_tag != tags[i];
}

void TfbsDPM::gibbsSample(unsigned int steps) {
        for (unsigned int i = 0; i < steps; i++) {
                double sum = 0;
                for (Data::iterator_randomized it = da->begin_randomized();
                     it != da->end_randomized(); it++) {
                        bool switched = sample(**it);
                        if (switched) sum+=1;
                }
                hist_switches.push_back(sum/da->size());
                compute_statistics();
        }
}
