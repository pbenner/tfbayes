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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "dpm.hh"

using namespace std;

DPM::DPM(Data* data)
        : da(data), cl(*da),
          // strength parameter for the dirichlet process
          alpha(1.0),
          // mixture weight for the dirichlet process
          lambda(0.1),
          tfbs_alpha(gsl_matrix_alloc(DPM::TFBS_LENGTH, DPM::NUCLEOTIDES)),
          bg_alpha(gsl_matrix_alloc(DPM::BG_LENGTH, DPM::NUCLEOTIDES))
{
        // initialize prior for the tfbs
        for (int i = 0; i < DPM::TFBS_LENGTH; i++) {
                for (int j = 0; j < DPM::NUCLEOTIDES; j++) {
                        gsl_matrix_set(tfbs_alpha, i, j, 1);
                }
        }
        // initialize prior for the background model
        for (int i = 0; i < DPM::BG_LENGTH; i++) {
                for (int j = 0; j < DPM::NUCLEOTIDES; j++) {
                        gsl_matrix_set(bg_alpha, i, j, DPM::TFBS_LENGTH);
                }
        }

        // initialize posterior predictive distributions for each cluster
        for (Cluster::iterator_all it = cl.begin_all(); it != cl.end_all(); it++) {
                if ((*it).tag == DPM::BG_CLUSTER) {
                        // background model
                        gsl_matrix* counts = gsl_matrix_alloc(DPM::BG_LENGTH, DPM::NUCLEOTIDES);
                        count_statistic(*it, bg_alpha, counts);
                        (*it).dist = new ProductDirichlet(counts);
                }
                else {
                        // tfbs models
                        gsl_matrix* counts = gsl_matrix_alloc(DPM::TFBS_LENGTH, DPM::NUCLEOTIDES);
                        count_statistic(*it, tfbs_alpha, counts);
                        (*it).dist = new ProductDirichlet(counts);
                }
        }

        // for sampling statistics
        hist_switches.push_back(0);
        hist_likelihood.push_back(likelihood());
}

DPM::~DPM() {
        delete(da);

        for (Cluster::iterator_all it = cl.begin_all(); it != cl.end_all(); it++) {
                delete((*it).dist);
        }
}

void
DPM::count_statistic(const Cluster::cluster& cluster, gsl_matrix* alpha, gsl_matrix* counts) {
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
                da->get_nucleotide(**it, len, buf);
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

double
DPM::likelihood() {
        return 0.0;
}

void
DPM::compute_statistics() {
}

bool
DPM::check_element(Data::element& element)
{
        // check if there is enough space for the binding site
        if (da->num_successors(element) < DPM::TFBS_LENGTH-1) {
                return false;
        }
        // check if this is alreade a binding site
        if (cl.getClusterTag(element) > 0) {
                return true;
        }
        if (cl.getClusterTag(element) == 0) {
                Data::iterator it = da->find(element);
                // check if all successing nucleotides belong to the background
                for (int i = 0; i < DPM::TFBS_LENGTH-1 && it != da->end(); i++) {
                        it++;
                        if (cl.getClusterTag(*it) != 0) {
                                return false;
                        }
                }
                // check if all previous nucleotides belong to the background
                it = da->find(element);
                for (int i = 0; i < DPM::TFBS_LENGTH-1 && it != da->begin(); i++) {
                        it--;
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
DPM::release_block(char* nucleotides, Data::element& element, Cluster::cluster& c) {
        // release a block of nucleotides from its clusters
        Data::iterator it = da->find(element);
        for (int i = 0; i < DPM::TFBS_LENGTH; i++) {
                cl.release(*it);
                it++;
        }
        if (c.tag == DPM::BG_CLUSTER) {
                for (int i = 0; i < DPM::TFBS_LENGTH; i++) {
                        ((ProductDirichlet*)c.dist)->remove_from_count_statistic(nucleotides+i);
                }
        }
        else {
                ((ProductDirichlet*)c.dist)->remove_from_count_statistic(nucleotides);
        }
}

void
DPM::assign_block(char* nucleotides, Data::element& element, Cluster::cluster& c) {
        if (c.tag == DPM::BG_CLUSTER) {
                // assign all nucleotides to the background cluster
                Data::iterator it = da->find(element);
                for (int i = 0; i < DPM::TFBS_LENGTH; i++) {
                        cl.assign(*it, DPM::BG_CLUSTER);
                        ((ProductDirichlet*)c.dist)->add_to_count_statistic(nucleotides+i);
                        it++;
                }
        }
        else {
                // this is a binding site:
                // assign this element to its class and leave
                // all remaining nucleotides unassigned
                cl.assign(element, c.tag);
                ((ProductDirichlet*)c.dist)->add_to_count_statistic(nucleotides);
        }
}

bool
DPM::sample(Data::element& element) {
        ////////////////////////////////////////////////////////////////////////
        // check if we can sample this element
        if (check_element(element) == false) {
                return false;
        }

        ////////////////////////////////////////////////////////////////////////
        // buffer to store the sequence of nucleotides, starting at `element'
        char nucleotides[DPM::TFBS_LENGTH];
        da->get_nucleotide(element, DPM::TFBS_LENGTH, nucleotides);

        ////////////////////////////////////////////////////////////////////////
        // release the element from its cluster
        Cluster::cluster_tag_t old_cluster_tag = cl.getClusterTag(element);
        release_block(nucleotides, element, cl[old_cluster_tag]);
        Cluster::size_type num_clusters = cl.size();
        double weights[num_clusters+1];
        Cluster::cluster_tag_t tags[num_clusters+1];
        double sum = 0;

        Cluster::cluster_tag_t i = 0;
        for (Cluster::iterator it = cl.begin(); it != cl.end(); it++) {
                tags[i] = (*it)->tag;
                ////////////////////////////////////////////////////////////////
                // mixture component 1: background model
                if (tags[i] == DPM::BG_CLUSTER) {
                        weights[i] = 1;
                        for (int j = 0; j < DPM::TFBS_LENGTH; j++) {
                                weights[i] *= (*it)->dist->pdf(nucleotides+j);
                        }
                        weights[i] *= (1-lambda);
                        sum += weights[i];
                }
                ////////////////////////////////////////////////////////////////
                // mixture component 2: dirichlet process for tfbs models
                else {
                        double num_elements = (double)(*it)->elements.size();
                        weights[i] = lambda*num_elements*(*it)->dist->pdf(nucleotides);
                        // normalization constant
                        sum += weights[i];
                }
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        tags[num_clusters]    = cl.next_free_cluster()->tag;
        weights[num_clusters] = alpha*cl[tags[num_clusters]].dist->pdf(nucleotides);
        sum += weights[num_clusters];

        ////////////////////////////////////////////////////////////////////////
        // normalize
        for (Cluster::cluster_tag_t i = 0; i < (Cluster::cluster_tag_t)num_clusters+1; i++) {
                weights[i] /= sum;
        }

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        gsl_ran_discrete_t* gdd  = gsl_ran_discrete_preproc(num_clusters+1, weights);
        i = gsl_ran_discrete(_r, gdd);
        gsl_ran_discrete_free(gdd);
        assign_block(nucleotides, element, cl[tags[i]]);

        ////////////////////////////////////////////////////////////////////////
        // return true if the cluster assignment has changed
        return old_cluster_tag != tags[i];
}

void
DPM::gibbsSample(unsigned int n) {
        // sample `n' times
        for (unsigned int i = 0; i < n; i++) {
                // loop through all elements
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

ostream& operator<< (ostream& o, DPM const& dpm)
{
        return o << dpm.cl;
}
