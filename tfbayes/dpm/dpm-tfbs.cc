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

#include <sstream>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <dpm-tfbs.hh>
#include <statistics.hh>

#include <tfbayes/logarithmetic.h>

#include <parsmm/abstract_set.h>
#include <parsmm/static_pars_tree.h>

using namespace std;

DpmTfbs::DpmTfbs(const tfbs_options_t& options, const data_tfbs_t& data)
        : // baseline
          _baseline_weights(options.baseline_weights),
          // raw sequences
          _data(data),
          // cluster manager
          _state(data.sizes(), options.tfbs_length, 0, _data),
          // mixture weight for the dirichlet process
          _lambda(options.lambda),
          _lambda_log(log(options.lambda)),
          _lambda_inv_log(log(1-options.lambda)),
          // length of tfbs
          _tfbs_length(options.tfbs_length)
{
        ////////////////////////////////////////////////////////////////////////////////
        // initialize samples
        for (size_t i = 0; i < data.size(); i++) {
                _samples.probabilities.push_back(vector<double>(data.size(i), 0.0));
        }

        ////////////////////////////////////////////////////////////////////////////////
        // add background model to the state
        if (options.background_model == "independence" || options.background_model == "") {
                const matrix<double>& bg_alpha = init_alpha(BG_LENGTH);
                ProductDirichlet* bg = new ProductDirichlet(bg_alpha, _data);
                bg_cluster_tag = _state.add_cluster(bg);
        }
        else if (options.background_model == "markov chain mixture") {
                assert(options.context >= 0);
                MarkovChainMixture* bg = new MarkovChainMixture(ALPHABET_SIZE, options.context, _data, _state.cluster_assignments, 0);
                bg_cluster_tag = _state.add_cluster(bg);
        }
        else if (options.background_model == "parsimonious tree") {
                assert(options.context >= 0);
                ParsimoniousTree* bg = new ParsimoniousTree(ALPHABET_SIZE, options.context, _data, _state.cluster_assignments, 0);
                bg_cluster_tag = _state.add_cluster(bg);
        }
        else {
                cerr << "Unknown background model." << endl;
                exit(EXIT_FAILURE);
        }
        // add model components for the baseline measure
        for (size_t i = 0; i < options.baseline_priors.size(); i++) {
                assert(options.tfbs_length == options.baseline_priors[i].size());
                ProductDirichlet* dirichlet = new ProductDirichlet(options.baseline_priors[i], _data);
                const model_tag_t model_tag = _state.add_baseline_model(dirichlet);
                _model_tags.push_back(model_tag);
        }

        ////////////////////////////////////////////////////////////////////////////////
        // assign all elements to the background
        for (da_iterator it = _data.begin();
             it != _data.end(); it++) {
                range_t range(**it, 1);
                _state[bg_cluster_tag].add_observations(range);
        }

        ////////////////////////////////////////////////////////////////////////////////
        // set the process prior
        if (options.process_prior == "pitman-yor process" || options.process_prior == "") {
                _process_prior = new pitman_yor_prior(_state, options.alpha, options.discount, bg_cluster_tag);
        }
        else if (options.process_prior == "uniform process") {
                _process_prior = new uniform_prior(_state, options.alpha);
        }
        else if (options.process_prior == "poppe process") {
                _process_prior = new poppe_prior(_state);
        }
        else {
                cerr << "Unknown prior process." << endl;
                exit(EXIT_FAILURE);
        }

        ////////////////////////////////////////////////////////////////////////////////
        //test();
        //test_background();
        //test_moves();
        //test_metropolis_hastings();
}

DpmTfbs::DpmTfbs(const DpmTfbs& dpm)
        : // baseline
          _baseline_weights(dpm._baseline_weights),
          _model_tags(dpm._model_tags),
          // raw sequences
          _data(dpm._data),
          // cluster manager
          _state(_state),
          // mixture weight for the dirichlet process
          _lambda(dpm._lambda),
          _lambda_log(dpm._lambda_log),
          _lambda_inv_log(dpm._lambda_inv_log),
          // length of tfbs
          _tfbs_length(dpm._tfbs_length),
          // samples
          _samples(dpm._samples),
          // process prios
          _process_prior(dpm._process_prior->clone())
{
}

DpmTfbs::~DpmTfbs() {
        delete(_process_prior);
}

DpmTfbs*
DpmTfbs::clone() const {
        return new DpmTfbs(*this);
}

bool
DpmTfbs::valid_for_sampling(const index_i& index) const
{
        const size_t sequence = index[0];
        const size_t position = index[1];

        // check if there is a tfbs starting here, if not check
        // succeeding positions
        if (_state.tfbs_start_positions[index] == 0) {
                // check if this element belongs to a tfbs that starts
                // earlier in the sequence
                if (_state[index] != bg_cluster_tag) {
                        return false;
                }
                if (_state.cluster_assignments[seq_index_t(sequence, position)] == -1) {
                        return false;
                }
                // check if there is a tfbs starting within the word
                for (size_t i = 1; i < _tfbs_length; i++) {
                        if (_state.tfbs_start_positions[seq_index_t(sequence, position+i)] == 1) {
                                return false;
                        }
                        if (_state.cluster_assignments[seq_index_t(sequence, position+i)] == -1) {
                                return false;
                        }
                }
        }

        return true;
}

size_t
DpmTfbs::mixture_components() const
{
        return _state.size();
}

size_t
DpmTfbs::baseline_components() const
{
        return _model_tags.size();
}

void
DpmTfbs::mixture_weights(const index_i& index, double log_weights[], cluster_tag_t cluster_tags[])
{
        const range_t range(index, _tfbs_length);
        ssize_t mixture_n  = mixture_components();
        ssize_t baseline_n = baseline_components();
        double sum         = -HUGE_VAL;

        cluster_tag_t i = 0;
        ////////////////////////////////////////////////////////////////////////
        // loop through existing clusters
        for (cm_iterator it = _state.begin(); it != _state.end(); it++) {
                cluster_t& cluster = **it;
                cluster_tags[i] = cluster.cluster_tag();
                if (cluster.cluster_tag() == bg_cluster_tag) {
                        ////////////////////////////////////////////////////////
                        // mixture component 1: background model
                        sum = logadd(sum, _lambda_inv_log + cluster.model().log_predictive(range));
                }
                else {
                        ////////////////////////////////////////////////////////
                        // mixture component 2: dirichlet process
                        sum = logadd(sum, _lambda_log + _process_prior->predictive(cluster) + cluster.model().log_predictive(range));
                }
                log_weights[i] = sum;
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (i = 0; i < baseline_n; i++) {
                cluster_t& cluster = _state.get_free_cluster(_model_tags[i]);
                cluster_tags[mixture_n+i] = cluster.cluster_tag();
                sum = logadd(sum, _lambda_log + _process_prior->predictive(cluster) + log(_baseline_weights[i]) +
                             cluster.model().log_predictive(range));
                log_weights[mixture_n+i] = sum;
        }
}

double
DpmTfbs::likelihood() const {
        double result = 0;

        for (cm_iterator it = _state.begin();
             it != _state.end(); it++) {
                cluster_t& cluster = **it;
                result += cluster.model().log_likelihood();
        }

        assert(!isnan(result));

        return result;
}

double
DpmTfbs::posterior() const {
        double result = likelihood();

        result += _state.num_tfbs*_lambda_log;
        result += _state[bg_cluster_tag].size()*_lambda_inv_log;
        result += _process_prior->joint();

        assert(!isnan(result));

        return result;
}

void
DpmTfbs::update_graph(sequence_data_t<short> tfbs_start_positions)
{
        // loop through all clusters
        for (cm_iterator it = _state.begin(); it != _state.end(); it++) {
                const cluster_t& cluster = **it;
                if (cluster.cluster_tag() != bg_cluster_tag) {
                        // loop through cluster elements
                        for (cl_iterator is = cluster.begin(); is != cluster.end(); is++) {
                                cl_iterator iu = is; iu++;
                                while (iu != cluster.end()) {
                                        // record edge
                                        _samples.graph[edge_t(is->index, iu->index)]++;
                                        iu++;
                                }
                        }
                }
        }
}

void
DpmTfbs::update_samples(size_t sampling_steps) {
        if (sampling_steps % 100 == 0) {
                _samples.graph.cleanup(1);
        }
        for (da_iterator it = _data.begin();
             it != _data.end(); it++) {
                const index_i& index  = **it;
                const size_t sequence = index[0];
                const size_t position = index[1];
                if (_state[index] == bg_cluster_tag) {
                        const double tmp   = _samples.probabilities[sequence][position];
                        const double value = ((double)sampling_steps*tmp)/((double)sampling_steps+1.0);
                        _samples.probabilities[sequence][position] = value;
                }
                else {
                        const double tmp   = _samples.probabilities[sequence][position];
                        const double value = ((double)sampling_steps*tmp+1.0)/((double)sampling_steps+1.0);
                        _samples.probabilities[sequence][position] = value;
                }
        }
        update_graph(_state.tfbs_start_positions);
}

samples_t&
DpmTfbs::samples() {
        return _samples;
}

dpm_tfbs_state_t&
DpmTfbs::state() {
        return _state;
}

const dpm_tfbs_state_t&
DpmTfbs::state() const {
        return _state;
}

std::matrix<double>
DpmTfbs::init_alpha(size_t length)
{
        std::matrix<double> alpha;

        // initialize prior for the background model
        for (size_t i = 0; i < length; i++) {
                alpha.push_back(std::vector<double>(ALPHABET_SIZE, 1));
        }
        return alpha;
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const DpmTfbs& dpm)
{
        o << "Cluster Assignments:"     << endl;
        o << dpm._state.cluster_assignments   << endl;

        o << "TFBS Start Positions:"    << endl;
        o << dpm._state.tfbs_start_positions  << endl;

        o << "Clusters:"                << endl;
        o << dpm._state;

        return o;
}
