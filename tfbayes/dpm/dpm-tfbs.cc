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

#define __process_prior ((*this).*(_process_prior))

DpmTfbs::DpmTfbs(const tfbs_options_t& options, const data_tfbs_t& data)
        : // length of tfbs
          TFBS_LENGTH(options.tfbs_length),
          // baseline
          _baseline_weights(options.baseline_weights),
          // raw sequences
          _data(data),
          // cluster manager
          _state(data.sizes(), TFBS_LENGTH, 0, _data),
          // strength parameter for the dirichlet process
          alpha(options.alpha),
          alpha_log(log(options.alpha)),
          // pitman-yor discount factor
          discount(options.discount),
          discount_log(log(options.discount)),
          // mixture weight for the dirichlet process
          lambda(options.lambda),
          lambda_log(log(options.lambda)),
          lambda_inv_log(log(1-options.lambda))
{
        ////////////////////////////////////////////////////////////////////////////////
        // initialize joint posterior
        for (size_t i = 0; i < data.size(); i++) {
                _posterior.probabilities.push_back(vector<double>(data.size(i), 0.0));
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
                _process_prior = &DpmTfbs::py_prior;
        }
        else if (options.process_prior == "uniform process") {
                _process_prior = &DpmTfbs::uniform_prior;
        }
        else if (options.process_prior == "poppe process") {
                _process_prior = &DpmTfbs::poppe_prior;
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

DpmTfbs::~DpmTfbs() {
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
                for (size_t i = 1; i < TFBS_LENGTH; i++) {
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
        const range_t range(index, TFBS_LENGTH);
        ssize_t mixture_n  = mixture_components();
        ssize_t baseline_n = baseline_components();
        double sum         = -HUGE_VAL;

        cluster_tag_t i = 0;
        ////////////////////////////////////////////////////////////////////////
        // loop through existing clusters
        for (cm_iterator it = _state.begin(); it != _state.end(); it++) {
                Cluster& cluster = **it;
                cluster_tags[i] = cluster.cluster_tag();
                if (cluster.cluster_tag() == bg_cluster_tag) {
                        ////////////////////////////////////////////////////////
                        // mixture component 1: background model
                        sum = logadd(sum, lambda_inv_log + cluster.model().log_predictive(range));
                }
                else {
                        ////////////////////////////////////////////////////////
                        // mixture component 2: dirichlet process
                        sum = logadd(sum, lambda_log + __process_prior(cluster) + cluster.model().log_predictive(range));
                }
                log_weights[i] = sum;
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (i = 0; i < baseline_n; i++) {
                Cluster& cluster = _state.get_free_cluster(_model_tags[i]);
                cluster_tags[mixture_n+i] = cluster.cluster_tag();
                sum = logadd(sum, lambda_log + __process_prior(cluster) + log(_baseline_weights[i]) +
                             cluster.model().log_predictive(range));
                log_weights[mixture_n+i] = sum;
        }
}

double
DpmTfbs::likelihood() const {
        double result = 0;

        for (cm_iterator it = _state.begin();
             it != _state.end(); it++) {
                Cluster& cluster = **it;
                result += cluster.model().log_likelihood();
        }
        return result;
}

bool
DpmTfbs::proposal(Cluster& cluster)
{
        double likelihood_ref = likelihood();
        double likelihood_new;

        _state.save(cluster.cluster_tag(), bg_cluster_tag);
        // propose move to right
        if (_state.move_right(cluster)) {
                likelihood_new = likelihood();

                if (likelihood_new > likelihood_ref) {
                        cout << "cluster "
                             << cluster.cluster_tag()
                             << ": move to right accepted"
                             << endl;
                        return true;
                }
        }

        // if not accepted, restore state
        _state.restore();
        // propose move to right
        if (_state.move_left(cluster)) {
                likelihood_new = likelihood();

                if (likelihood_new > likelihood_ref) {
                        cout << "cluster "
                             << cluster.cluster_tag()
                             << ": move to left accepted"
                             << endl;
                        return true;
                }
        }

        // if not accepted, restore state
        _state.restore();

        return false;
}

void
DpmTfbs::metropolis_hastings()
{
        for (cm_iterator it = _state.begin(); it != _state.end(); it++) {
                Cluster& cluster = **it;
                if (cluster.cluster_tag() != bg_cluster_tag && cluster.size() > 1) {
                        proposal(cluster);
                }
        }
}

void
DpmTfbs::update_graph(sequence_data_t<short> tfbs_start_positions)
{
        // loop through all clusters
        for (cm_iterator it = _state.begin(); it != _state.end(); it++) {
                const Cluster& cluster = **it;
                if (cluster.cluster_tag() != bg_cluster_tag) {
                        // loop through cluster elements
                        for (cl_iterator is = cluster.begin(); is != cluster.end(); is++) {
                                cl_iterator iu = is; iu++;
                                while (iu != cluster.end()) {
                                        // record edge
                                        _posterior.graph[edge_t(is->index, iu->index)]++;
                                        iu++;
                                }
                        }
                }
        }
}

void
DpmTfbs::update_posterior(size_t sampling_steps) {
        metropolis_hastings();
        if (sampling_steps % 100 == 0) {
                _posterior.graph.cleanup(1);
        }
        for (da_iterator it = _data.begin();
             it != _data.end(); it++) {
                const index_i& index  = **it;
                const size_t sequence = index[0];
                const size_t position = index[1];
                if (_state[index] == bg_cluster_tag) {
                        const double tmp   = _posterior.probabilities[sequence][position];
                        const double value = ((double)sampling_steps*tmp)/((double)sampling_steps+1.0);
                        _posterior.probabilities[sequence][position] = value;
                }
                else {
                        const double tmp   = _posterior.probabilities[sequence][position];
                        const double value = ((double)sampling_steps*tmp+1.0)/((double)sampling_steps+1.0);
                        _posterior.probabilities[sequence][position] = value;
                }
        }
        update_graph(_state.tfbs_start_positions);
}

posterior_t&
DpmTfbs::posterior() {
        return _posterior;
}

const data_tfbs_t&
DpmTfbs::data() const {
        return _data;
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
