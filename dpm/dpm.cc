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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <dpm.hh>

using namespace std;

DPM::DPM(const Data& data)
        : _data(data),
          // strength parameter for the dirichlet process
          alpha(0.07),
          // mixture weight for the dirichlet process
          lambda(0.02),
          // priors
          bg_alpha(gsl_matrix_alloc(DPM::BG_LENGTH, DPM::NUCLEOTIDES)),
          tfbs_alpha(gsl_matrix_alloc(DPM::TFBS_LENGTH, DPM::NUCLEOTIDES)),
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
        for (size_t i = 0; i < data.length(); i++) {
                _posterior.push_back(vector<double>(data.length(i), 0.0));
        }

        // starting positions of tfbs
        for(size_t i = 0; i < data.length(); i++) {
                this->tfbs_start_positions.push_back(vector<bool>(data.length(i), false));
        }

        // initialize cluster manager
        ProductDirichlet* tfbs_product_dirichlet = new ProductDirichlet(tfbs_alpha);
        ProductDirichlet* bg_product_dirichlet   = new ProductDirichlet(bg_alpha);
        _cluster_manager = new ClusterManager(_data, tfbs_product_dirichlet);
        bg_cluster_tag   = _cluster_manager->add_cluster(bg_product_dirichlet);

        // assign all elements to the background
        for (Data::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                const word_t word = _data.get_word(*it, DPM::BG_LENGTH);
                (*_cluster_manager)[bg_cluster_tag].add_word(word);
        }
}

DPM::~DPM() {
        gsl_matrix_free(tfbs_alpha);
        gsl_matrix_free(bg_alpha);

        delete(_cluster_manager);
}

bool
DPM::valid_for_sampling(const element_t& element, const word_t& word)
{
        const size_t sequence = word.sequence;
        const size_t position = word.position;
        const size_t length   = word.length;

        // check if there is enough space
        if (_data.length(sequence) - position < length) {
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

double
DPM::compute_likelihood() {
        return 0.0;
}

void
DPM::update_posterior(size_t sampling_steps) {
        for (Data::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                const element_t& element = *it;
                const size_t sequence    = element.sequence;
                const size_t position    = element.position;
                if (_cluster_manager->get_cluster_tag(element) == bg_cluster_tag) {
                        double tmp   = _posterior[sequence][position];
                        double value = (sampling_steps*tmp)/(sampling_steps+1.0);
                        _posterior[sequence][position] = value;
                }
                else {
                        double tmp   = _posterior[sequence][position];
                        double value = (sampling_steps*tmp+1.0)/(sampling_steps+1.0);
                        _posterior[sequence][position] = value;
                }
        }
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const DPM& dpm)
{
        for (size_t i = 0; i < dpm._data.length(); i++) {
                for (size_t j = 0; j < dpm._data.length(i); j++) {
                        o << dpm.tfbs_start_positions[i][j] << " ";
                }
                o << endl;
        }
        return o;
}
