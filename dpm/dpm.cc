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
                _posterior.push_back(vector<double>(length, 0.0));
        }

        // initialize data structures for the gibbs sampler
        for(size_t i = 0; i < n; i++) {
                size_t m = strlen(sequences[i]);
                this->tfbs_start_positions.push_back(vector<bool>(m, false));
        }

        // initialize cluster manager
        ProductDirichlet* tfbs_product_dirichlet = new ProductDirichlet(tfbs_alpha);
        ProductDirichlet* bg_product_dirichlet   = new ProductDirichlet(bg_alpha);
        _data            = new Data(n, sequences);
        _cluster_manager = new ClusterManager(*_data, tfbs_product_dirichlet);
        bg_cluster_tag   = _cluster_manager->add_cluster(bg_product_dirichlet);

        // assign all elements to the background
        for (Data::iterator it = _data->begin();
             it != _data->end(); it++) {
                const word_t word = _data->get_word(*it, DPM::BG_LENGTH);
                (*_cluster_manager)[bg_cluster_tag].add_word(word);
        }

        // for sampling statistics
        hist_switches.push_back(0);
        hist_likelihood.push_back(compute_likelihood());
}

DPM::~DPM() {
        gsl_matrix_free(tfbs_alpha);
        gsl_matrix_free(bg_alpha);

        delete(_data);
        delete(_cluster_manager);
}

bool
DPM::valid_for_sampling(const element_t& element, const word_t& word)
{
        const size_t sequence = word.sequence;
        const size_t position = word.position;
        const size_t length   = word.length;

        // check if there is enough space
        if (_data->length(sequence) - position < length) {
                return false;
        }
        // check if there is a tfbs starting here, if not check
        // succeeding positions
        if (tfbs_start_positions[element.sequence][element.position] == 0) {
                // check if this element belongs to a tfbs that starts
                // earlier in the sequence
                if ((*_cluster_manager)[element] != bg_cluster_tag) {
                        return false;
                }
                // check if there is a tfbs starting within the word
                for (size_t i = 1; i < length; i++) {
                        if (tfbs_start_positions[sequence][position+i] == 1) {
                                return false;
                        }
                }
        }

        return true;
}

void
DPM::add_word(const word_t& word, cluster_tag_t tag)
{
        (*_cluster_manager)[tag].add_word(word);
        if (tag != bg_cluster_tag) {
                num_tfbs++;
                tfbs_start_positions[word.sequence][word.position] = 1;
        }
}

void
DPM::remove_word(const word_t& word, cluster_tag_t tag)
{
        (*_cluster_manager)[tag].remove_word(word);
        if (tag != bg_cluster_tag) {
                num_tfbs--;
                tfbs_start_positions[word.sequence][word.position] = 0;
        }
}

size_t
DPM::mixture_components() const
{
        return _cluster_manager->size() + 1;
}

void
DPM::mixture_weights(const word_t& word, double weights[], cluster_tag_t tags[]) const
{
        size_t components = mixture_components();
        double dp_norm    = num_tfbs + alpha;
        double sum        = 0;

        cluster_tag_t i = 0;
        for (ClusterManager::iterator it = _cluster_manager->begin(); it != _cluster_manager->end(); it++) {
                Cluster& cluster = **it;
                tags[i] = cluster.tag();
                ////////////////////////////////////////////////////////////////
                // mixture component 1: background model
                if (tags[i] == bg_cluster_tag) {
                        weights[i] = (1-lambda)*cluster.distribution().pdf(word);
                        // normalization constant
                        sum += weights[i];
                }
                ////////////////////////////////////////////////////////////////
                // mixture component 2: dirichlet process for tfbs models
                else {
                        double num_elements = (double)cluster.size();
                        weights[i] = lambda*num_elements/dp_norm*cluster.distribution().pdf(word);
                        // normalization constant
                        sum += weights[i];
                }
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        tags[components-1] = _cluster_manager->get_free_cluster().tag();
        weights[components-1] = alpha/dp_norm*(*_cluster_manager)[tags[components-1]].distribution().pdf(word);
        sum += weights[components-1];

        ////////////////////////////////////////////////////////////////////////
        // normalize
        for (size_t i = 0; i < components; i++) {
                weights[i] /= sum;
        }
}

bool
DPM::sample(const element_t& element) {
        const word_t word = _data->get_word(element, DPM::TFBS_LENGTH);
        ////////////////////////////////////////////////////////////////////////
        // check if we can sample this element
        if (!valid_for_sampling(element, word)) {
                return false;
        }
        ////////////////////////////////////////////////////////////////////////
        // release the element from its cluster
        cluster_tag_t old_cluster_tag = _cluster_manager->get_cluster_tag(element);
        remove_word(word, old_cluster_tag);
        size_t components = mixture_components();
        double weights[components];
        cluster_tag_t tags[components];
        mixture_weights(word, weights, tags);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        gsl_ran_discrete_t* gdd  = gsl_ran_discrete_preproc(components, weights);
        cluster_tag_t i = gsl_ran_discrete(_r, gdd);
        gsl_ran_discrete_free(gdd);
        cluster_tag_t new_cluster_tag = tags[i];

        ////////////////////////////////////////////////////////////////////////
        add_word(word, new_cluster_tag);

        return old_cluster_tag != new_cluster_tag;
}

// sampling methods
////////////////////////////////////////////////////////////////////////////////

double
DPM::compute_likelihood() {
        return 0.0;
}

void
DPM::update_posterior() {
        for (Data::iterator it = _data->begin();
             it != _data->end(); it++) {
                const element_t& element = *it;
                const size_t sequence    = element.sequence;
                const size_t position    = element.position;
                if (_cluster_manager->get_cluster_tag(element) == bg_cluster_tag) {
                        double tmp   = _posterior[sequence][position];
                        double value = (total_sampling_steps*tmp)/(total_sampling_steps+1.0);
                        _posterior[sequence][position] = value;
                }
                else {
                        double tmp   = _posterior[sequence][position];
                        double value = (total_sampling_steps*tmp+1.0)/(total_sampling_steps+1.0);
                        _posterior[sequence][position] = value;
                }
        }
}

void
DPM::gibbs_sample(size_t n, size_t burnin) {
        // burn in sampling
        for (size_t i = 0; i < burnin; i++) {
                printf("Burn in... [%u]\n", (unsigned int)i+1);
                for (Data::iterator_randomized it = _data->begin_randomized();
                     it != _data->end_randomized(); it++) {
                        sample(**it);
                }
        }
        // sample `n' times
        for (size_t i = 0; i < n; i++) {
                // loop through all elements
                printf("Sampling... [%u]\n", (unsigned int)i+1);
                double sum = 0;
                for (Data::iterator_randomized it = _data->begin_randomized();
                     it != _data->end_randomized(); it++) {
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
        for (size_t i = 0; i < dpm._data->length(); i++) {
                for (size_t j = 0; j < dpm._data->length(i); j++) {
                        o << dpm.tfbs_start_positions[i][j] << " ";
                }
                o << endl;
        }
        return o;
}
