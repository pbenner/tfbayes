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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <dpm.hh>

using namespace std;

DPM::DPM(size_t n, char *sequences[])
        : // strength parameter for the dirichlet process
          alpha(0.07),
          // mixture weight for the dirichlet process
          lambda(0.02),
          // priors
          bg_alpha(gsl_matrix_alloc(DPM::BG_LENGTH, DPM::NUCLEOTIDES)),
          tfbs_alpha(gsl_matrix_alloc(DPM::TFBS_LENGTH, DPM::NUCLEOTIDES)),
          // sampling history and posterior distribution
          total_sampling_steps(0),
          // number of transcription factor binding sites
          num_tfbs(0)
{
        // initialize prior for the tfbs
        for (size_t i = 0; i < DPM::TFBS_LENGTH; i++) {
                for (size_t j = 0; j < DPM::NUCLEOTIDES; j++) {
                        gsl_matrix_set(tfbs_alpha, i, j, 1);
                }
        }

        // initialize prior for the background model
        for (size_t i = 0; i < DPM::BG_LENGTH; i++) {
                for (size_t j = 0; j < DPM::NUCLEOTIDES; j++) {
                        gsl_matrix_set(bg_alpha, i, j, 1);
                }
        }

        // initialize joint posterior
        for (size_t i = 0; i < n; i++) {
                const size_t length = strlen(sequences[i]);
                posterior.push_back(vector<double>(length, 0.0));
        }

        // initialize data structures for the gibbs sampler
        for(size_t i = 0; i < n; i++) {
                size_t m = strlen(sequences[i]);
                this->blocked_for_sampling.push_back(vector<bool>(m, false));
        }

        // initialize cluster manager
        ProductDirichlet* tfbs_product_dirichlet = new ProductDirichlet(tfbs_alpha);
        ProductDirichlet* bg_product_dirichlet   = new ProductDirichlet(bg_alpha);
        cluster_manager = new ClusterManager(tfbs_product_dirichlet);
        bg_cluster_tag  = cluster_manager->add_cluster(bg_product_dirichlet);
        data            = new Data(n, sequences, bg_cluster_tag);

        // for sampling statistics
        hist_switches.push_back(0);
        hist_likelihood.push_back(compute_likelihood());
}

DPM::~DPM() {
        gsl_matrix_free(tfbs_alpha);
        gsl_matrix_free(bg_alpha);

        delete(data);
        delete(cluster_manager);
}

void
DPM::block_for_sampling(const word_t& word)
{
        const size_t sequence = word.sequence;
        const size_t position = word.position;
        const size_t length   = word.length;

        for (size_t i = 1; i < length; i++) {
                blocked_for_sampling[sequence][position+i] = true;
        }
}

void
DPM::unblock_for_sampling(const word_t& word)
{
        const size_t sequence = word.sequence;
        const size_t position = word.position;
        const size_t length   = word.length;

        for (size_t i = 1; i < length; i++) {
                blocked_for_sampling[sequence][position+i] = false;
        }
}

bool
DPM::valid_for_sampling(const word_t& word)
{
        const size_t sequence = word.sequence;
        const size_t position = word.position;
        const size_t length   = word.length;

        // check if there is enough space
        if (data->length(sequence) - position < length) {
                return false;
        }
        // check if a position within the word is blocked
        for (size_t i = 0; i < length; i++) {
                if (blocked_for_sampling[sequence][position+i]) {
                        return false;
                }
        }
        if (data->length(sequence) - position < length-1) {
                if (blocked_for_sampling[sequence][position+length]) {
                        return false;
                }
        }

        return true;
}

bool
DPM::sample(const element_t& element) {
        word_t word;
        data->get_word(element, DPM::TFBS_LENGTH, word);
        ////////////////////////////////////////////////////////////////////////
        // check if we can sample this element
        if (!valid_for_sampling(word)) {
                return false;
        }
        cout << word << endl;
        cout << *data << endl;
        ////////////////////////////////////////////////////////////////////////
        // release the element from its cluster
        cluster_tag_t old_cluster_tag = data->get_cluster_tag(element);
        if (old_cluster_tag == bg_cluster_tag) {
                printf("Sampling background\n");
        }
        else {
                printf("Sampling tfbs\n");
        }
        (*cluster_manager)[old_cluster_tag].remove_word(word);
        size_t num_clusters = cluster_manager->size();
        printf("num_tfbs: %d\n", (int)num_tfbs);
        printf("num_clusters: %d\n", (int)num_clusters);
        double weights[num_clusters+1];
        cluster_tag_t tags[num_clusters+1];
        double dp_norm = num_tfbs + alpha;
        double sum = 0;

        cluster_tag_t i = 0;
        for (ClusterManager::iterator it = cluster_manager->begin(); it != cluster_manager->end(); it++) {
                Cluster& cluster = **it;
                tags[i] = cluster.tag;
                ////////////////////////////////////////////////////////////////
                // mixture component 1: background model
                if (tags[i] == bg_cluster_tag) {
                        weights[i] = (1-lambda)*cluster.distribution->pdf(word);
                        // normalization constant
                        sum += weights[i];
                }
                ////////////////////////////////////////////////////////////////
                // mixture component 2: dirichlet process for tfbs models
                else {
                        double num_elements = (double)cluster.size();
                        printf("num_elements: %d\n", (int)num_elements);
                        weights[i] = lambda*num_elements/dp_norm*cluster.distribution->pdf(word);
                        // normalization constant
                        sum += weights[i];
                }
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        tags[num_clusters]    = cluster_manager->get_free_cluster().tag;
        printf("free cluster tag: %d\n", (int)tags[num_clusters]);
        weights[num_clusters] = alpha/dp_norm*(*cluster_manager)[tags[num_clusters]].distribution->pdf(word);
        sum += weights[num_clusters];

        ////////////////////////////////////////////////////////////////////////
        // normalize
        for (size_t i = 0; i < num_clusters+1; i++) {
                if (tags[i] == bg_cluster_tag) {
                        printf("     weight[%d](bg): %f\n", (int)i, (float)weights[i]);
                }
                else {
                        printf("     weight[%d](tfbs): %f\n", (int)i, (float)weights[i]);
                }
                weights[i] /= sum;
                if (tags[i] == bg_cluster_tag) {
                        printf("NORM weight[%d](bg): %f\n", (int)i, (float)weights[i]);
                }
                else {
                        printf("NORM weight[%d](tfbs): %f\n", (int)i, (float)weights[i]);
                }
        }

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        gsl_ran_discrete_t* gdd  = gsl_ran_discrete_preproc(num_clusters+1, weights);
        i = gsl_ran_discrete(_r, gdd);
        gsl_ran_discrete_free(gdd);
        cluster_tag_t new_cluster_tag = tags[i];

        ////////////////////////////////////////////////////////////////////////
        // if the cluster assignment has changes, record it and return true
        if (old_cluster_tag == new_cluster_tag) {
                (*cluster_manager)[old_cluster_tag].add_word(word);
                cout << "Keeping cluster assignment" << endl;
                cout << endl << endl;
                return false;
        }
        else {
                if (old_cluster_tag == bg_cluster_tag && new_cluster_tag != bg_cluster_tag) {
                        block_for_sampling(word);
                        num_tfbs++;
                        cout << "BG -> TFBS" << endl;
                        cout << endl << endl;
                }
                else if (old_cluster_tag != bg_cluster_tag && new_cluster_tag == bg_cluster_tag) {
                        unblock_for_sampling(word);
                        num_tfbs--;
                        cout << "TFBS -> BG" << endl;
                        cout << endl << endl;
                }
                else {
                        cout << "TFBS -> TFBS" << endl;
                        cout << endl << endl;
                }
                (*cluster_manager)[new_cluster_tag].add_word(word);
                data->record_cluster_assignment(word, new_cluster_tag);
                return true;
        }
}

// sampling methods
////////////////////////////////////////////////////////////////////////////////

double
DPM::compute_likelihood() {
        return 0.0;
}

void
DPM::update_posterior() {
        for (Data::iterator it = data->begin();
             it != data->end(); it++) {
                const element_t& element = *it;
                const size_t sequence    = element.sequence;
                const size_t position    = element.position;
                if (data->get_cluster_tag(element) == bg_cluster_tag) {
                        double tmp   = posterior[sequence][position];
                        double value = (total_sampling_steps*tmp)/(total_sampling_steps+1.0);
                        posterior[sequence][position] = value;
                }
                else {
                        double tmp   = posterior[sequence][position];
                        double value = (total_sampling_steps*tmp+1.0)/(total_sampling_steps+1.0);
                        posterior[sequence][position] = value;
                }
        }
}

void
DPM::gibbs_sample(size_t n, size_t burnin) {
        // burn in sampling
        for (size_t i = 0; i < burnin; i++) {
                printf("Burn in... [%u]\n", (unsigned int)i+1);
                for (Data::iterator_randomized it = data->begin_randomized();
                     it != data->end_randomized(); it++) {
//                        printf("Burn in... [%u:%lu:%lu:%lu]\n", i+1, (*it)->x[0], (*it)->x[1], cl.size());
                        sample(**it);
                }
        }
        // sample `n' times
        for (size_t i = 0; i < n; i++) {
                // loop through all elements
                printf("Sampling... [%u]\n", (unsigned int)i+1);
                double sum = 0;
                for (Data::iterator_randomized it = data->begin_randomized();
                     it != data->end_randomized(); it++) {
                        bool switched = sample(**it);
                        if (switched) sum+=1;
                }
                hist_switches.push_back(sum);
                update_posterior();
                total_sampling_steps++;
        }
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const DPM& dpm)
{
        for (size_t i = 0; i < dpm.data->length(); i++) {
                for (size_t j = 0; j < dpm.data->length(i); j++) {
                        o << dpm.blocked_for_sampling[i][j] << " ";
                }
                o << endl;
        }
        return o;
}
