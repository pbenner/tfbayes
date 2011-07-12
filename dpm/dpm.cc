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
#include <tfbayes/logarithmetic.h>

using namespace std;

DPM::DPM(double alpha, double lambda, size_t tfbs_length, const Data& data)
        : // length of tfbs
          TFBS_LENGTH(tfbs_length),
          // priors
          bg_alpha(init_alpha(BG_LENGTH)),
          tfbs_alpha(init_alpha(TFBS_LENGTH)),
          // raw sequences
          _data(data),
          _cluster_manager(data, new ProductDirichlet(tfbs_alpha)),
          // strength parameter for the dirichlet process
          alpha(alpha),
          // mixture weight for the dirichlet process
          lambda(lambda),
          // number of transcription factor binding sites
          num_tfbs(0)
{
        // initialize joint posterior
        for (size_t i = 0; i < data.length(); i++) {
                _posterior.push_back(vector<double>(data.length(i), 0.0));
        }

        // starting positions of tfbs
        for(size_t i = 0; i < data.length(); i++) {
                this->tfbs_start_positions.push_back(vector<bool>(data.length(i), false));
        }

        // initialize cluster manager
        ProductDirichlet* bg_product_dirichlet   = new ProductDirichlet(bg_alpha);
        bg_cluster_tag   = _cluster_manager.add_cluster(bg_product_dirichlet);

        // assign all elements to the background
        for (Data::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                const word_t word = _data.get_word(*it, BG_LENGTH);
                _cluster_manager[bg_cluster_tag].add_word(word);
        }
}

DPM::~DPM() {
        gsl_matrix_free(tfbs_alpha);
        gsl_matrix_free(bg_alpha);
}

DPM*
DPM::clone() const {
        return new DPM(*this);
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
                if (_cluster_manager[element] != bg_cluster_tag) {
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
        _cluster_manager[tag].add_word(word);
        if (tag != bg_cluster_tag) {
                num_tfbs++;
                tfbs_start_positions[word.sequence][word.position] = 1;
        }
}

void
DPM::remove_word(const word_t& word, cluster_tag_t tag)
{
        _cluster_manager[tag].remove_word(word);
        if (tag != bg_cluster_tag) {
                num_tfbs--;
                tfbs_start_positions[word.sequence][word.position] = 0;
        }
}

size_t
DPM::mixture_components() const
{
        return _cluster_manager.size() ;
}

void
DPM::mixture_weights(const word_t& word, double weights[], cluster_tag_t tags[])
{
        size_t components = mixture_components();
        double dp_norm    = num_tfbs + alpha;
        double sum        = -HUGE_VAL;

        cluster_tag_t i = 0;
        for (ClusterManager::iterator it = _cluster_manager.begin(); it != _cluster_manager.end(); it++) {
                Cluster& cluster = **it;
                tags[i] = cluster.tag();
                ////////////////////////////////////////////////////////////////
                // mixture component 1: background model
                if (tags[i] == bg_cluster_tag) {
                        weights[i] = log(1-lambda) + cluster.distribution().log_pdf(word);
                        // normalization constant
                        sum = logadd(sum, weights[i]);
                }
                ////////////////////////////////////////////////////////////////
                // mixture component 2: dirichlet process for tfbs models
                else {
                        double num_elements = (double)cluster.size();
                        weights[i] = log(lambda*num_elements/dp_norm) + cluster.distribution().log_pdf(word);
                        // normalization constant
                        sum = logadd(sum, weights[i]);
                }
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        tags[components]    = _cluster_manager.get_free_cluster().tag();
        weights[components] = log(lambda*alpha/dp_norm) + _cluster_manager[tags[components]].distribution().log_pdf(word);
        sum = logadd(sum, weights[components]);

        ////////////////////////////////////////////////////////////////////////
        // normalize
        for (size_t i = 0; i < components+1; i++) {
                weights[i] = exp(weights[i] - sum);
        }
}

double
DPM::likelihood() const {
        double result = 0;

        for (ClusterManager::const_iterator it = _cluster_manager.begin();
             it != _cluster_manager.end(); it++) {
                Cluster& cluster = **it;
                result += cluster.distribution().log_likelihood();
        }
        return result;
}

void
DPM::update_posterior(size_t sampling_steps) {
        for (Data::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                const element_t& element = *it;
                const size_t sequence    = element.sequence;
                const size_t position    = element.position;
                if (_cluster_manager.get_cluster_tag(element) == bg_cluster_tag) {
                        double tmp   = _posterior[sequence][position];
                        double value = ((double)sampling_steps*tmp)/((double)sampling_steps+1.0);
                        _posterior[sequence][position] = value;
                }
                else {
                        double tmp   = _posterior[sequence][position];
                        double value = ((double)sampling_steps*tmp+1.0)/((double)sampling_steps+1.0);
                        _posterior[sequence][position] = value;
                }
        }
}

const posterior_t&
DPM::posterior() const {
        return _posterior;
}

const Data&
DPM::data() const {
        return _data;
}

const ClusterManager&
DPM::cluster_manager() const {
        return _cluster_manager;
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
