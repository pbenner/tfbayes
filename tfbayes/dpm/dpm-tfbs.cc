/* Copyright (C) 2011, 2012 Philipp Benner
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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <sstream>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/dpm/statistics.hh>
#include <tfbayes/utility/logarithmetic.h>

using namespace std;

dpm_tfbs_t::dpm_tfbs_t(const tfbs_options_t& options, const data_tfbs_t& data, const alignment_set_t<short>& alignment_set)
        : // baseline
          _baseline_weights(*options.baseline_weights),
          // phylogenetic information
          _data(data),
          // coded nucleotide sequences
          _alignment_set(alignment_set),
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
        // check that the alignment data matches the phylogenetic data
        assert(data.size() == alignment_set.size());
        for (size_t i = 0; i < data.size(); i++) {
                assert(data[i].size() == alignment_set[i].length());
        }

        ////////////////////////////////////////////////////////////////////////////////
        // initialize samples
        for (size_t i = 0; i < data.size(); i++) {
                _samples.probabilities.push_back(vector<double>(data.size(i), 0.0));
        }
        _samples.map_value = -HUGE_VAL;

        ////////////////////////////////////////////////////////////////////////////////
        // add background model to the state
        if (*options.background_model == "independence-dirichlet" || *options.background_model == "") {
                /* every position in the background is fully
                 * independet, this give more flexibility to the
                 * prior pseudocounts */
                independence_background_t* bg = new independence_background_t(*options.background_alpha, _data, _state.cluster_assignments());
                _state.bg_cluster_tag = _state.add_cluster(bg);
                bg->set_bg_cluster_tag(_state.bg_cluster_tag);
        }
        else if (*options.background_model == "independence-dirichlet-gamma") {
                /* every position in the background is fully
                 * independet, where the Dirichlet pseudocounts
                 * are integrated out */
                cerr << "TODO: variable gamma prior parameters" << endl;
                independence_background_t* bg = new independence_background_t(2.0, 2.0, _data, _state.cluster_assignments());
                _state.bg_cluster_tag = _state.add_cluster(bg);
                bg->set_bg_cluster_tag(_state.bg_cluster_tag);
        }
        else if (*options.background_model == "dirichlet") {
                /* single dirichlet-compound distribution for all
                 * nucleotides in the background */
                product_dirichlet_t* bg = new product_dirichlet_t(*options.background_alpha, _data);
                _state.bg_cluster_tag = _state.add_cluster(bg);
        }
        else if (*options.background_model == "uniform") {
                /* all sequences have the same probability (no phylogeny!) */
                uniform_background_t* bg = new uniform_background_t(_data, _alignment_set);
                _state.bg_cluster_tag = _state.add_cluster(bg);
        }
        else if (*options.background_model == "markov chain mixture") {
                assert(options.background_context >= 0);
                markov_chain_mixture_t* bg = new markov_chain_mixture_t(data_tfbs_t::alphabet_size, options, _data, _state.cluster_assignments(), 0);
                _state.bg_cluster_tag = _state.add_cluster(bg);
        }
        else {
                cerr << "Unknown background model." << endl;
                exit(EXIT_FAILURE);
        }
        // baseline weights are already initialized
        for (size_t i = 0; i < options.baseline_n; i++) {
                assert((*options.baseline_priors)[i].size() == options.tfbs_length);
                for (size_t j = 0; j < options.tfbs_length; j++) {
                        assert((*options.baseline_priors)[i][j].size() == data_tfbs_t::alphabet_size);
                }
                product_dirichlet_t* dirichlet = new product_dirichlet_t((*options.baseline_priors)[i], _data);
                _state.add_baseline_model(dirichlet, (*options.baseline_tags)[i]);
                _baseline_tags.push_back((*options.baseline_tags)[i]);
        }

        ////////////////////////////////////////////////////////////////////////////////
        // assign all elements to the background
        for (da_iterator it = _data.begin();
             it != _data.end(); it++) {
                // cannot use _state.add() here because we need to add
                // sites of length one to the background
                range_t range(**it, 1);
                _state[_state.bg_cluster_tag].add_observations(range);
        }

        ////////////////////////////////////////////////////////////////////////////////
        // set the process prior
        if (*options.process_prior == "pitman-yor process" || *options.process_prior == "") {
                _process_prior = new pitman_yor_prior(options.alpha, options.discount, _state.bg_cluster_tag);
        }
        else if (*options.process_prior == "uniform process") {
                _process_prior = new uniform_prior(options.alpha);
        }
        else if (*options.process_prior == "poppe process") {
                _process_prior = new poppe_prior();
        }
        else {
                cerr << "Unknown prior process." << endl;
                exit(EXIT_FAILURE);
        }

        ////////////////////////////////////////////////////////////////////////////////
        // use a map partition from a previous sampling run to
        // initialize the state
        if (options.partition) {
                for (dpm_partition_t::const_iterator it = options.partition->begin(); it != options.partition->end(); it++) {
                        const dpm_subset_t& subset(*it);
                        cluster_t& cluster = _state.get_free_cluster(subset.dpm_subset_tag());

                        for (dpm_subset_t::const_iterator is = subset.begin(); is != subset.end(); is++) {
                                assert(valid_for_sampling(**is));
                                assert(_state[**is] == _state.bg_cluster_tag);
                                _state.remove(**is, _state.bg_cluster_tag);
                                _state.add(**is, cluster.cluster_tag());
                        }
                }
        }
        ////////////////////////////////////////////////////////////////////////////////
        //test();
        //test_background();
        //test_moves();
        //test_metropolis_hastings();
}

dpm_tfbs_t::dpm_tfbs_t(const dpm_tfbs_t& dpm)
        : // baseline
          _baseline_weights(dpm._baseline_weights),
          _baseline_tags(dpm._baseline_tags),
          // phylogenetic data
          _data(dpm._data),
          // coded nucleotide sequences
          _alignment_set(dpm._alignment_set),
          // cluster manager
          _state(dpm._state),
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
{ }

dpm_tfbs_t::~dpm_tfbs_t() {
        delete(_process_prior);
}

dpm_tfbs_t*
dpm_tfbs_t::clone() const {
        return new dpm_tfbs_t(*this);
}

bool
dpm_tfbs_t::valid_for_sampling(const index_i& index) const
{
        const size_t sequence = index[0];
        const size_t position = index[1];

        // check if there is a tfbs starting here, if not check
        // succeeding positions
        if (_state.tfbs_start_positions[index] == 0) {
                // check if this element belongs to a tfbs that starts
                // earlier in the sequence
                if (_state[index] != _state.bg_cluster_tag) {
                        return false;
                }
                if (_state.cluster_assignments()[seq_index_t(sequence, position)] == -1) {
                        return false;
                }
                // check if there is a tfbs starting within the word
                for (size_t i = 1; i < _tfbs_length; i++) {
                        if (_state.tfbs_start_positions[seq_index_t(sequence, position+i)] == 1) {
                                return false;
                        }
                        if (_state.cluster_assignments()[seq_index_t(sequence, position+i)] == -1) {
                                return false;
                        }
                }
        }

        return true;
}

size_t
dpm_tfbs_t::mixture_components() const
{
        return _state.size();
}

size_t
dpm_tfbs_t::baseline_components() const
{
        return _baseline_tags.size();
}

void
dpm_tfbs_t::mixture_weights(const index_i& index, double log_weights[], cluster_tag_t cluster_tags[], const double temp)
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
                if (cluster.cluster_tag() == _state.bg_cluster_tag) {
                        ////////////////////////////////////////////////////////
                        // mixture component 1: background model
                        sum = logadd(sum, (_lambda_inv_log + cluster.model().log_predictive(range))/temp);
                }
                else {
                        ////////////////////////////////////////////////////////
                        // mixture component 2: dirichlet process
                        sum = logadd(sum, (_lambda_log + _process_prior->log_predictive(cluster, _state) + cluster.model().log_predictive(range))/temp);
                }
                log_weights[i] = sum;
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (i = 0; i < baseline_n; i++) {
                cluster_t& cluster = _state.get_free_cluster(_baseline_tags[i]);
                cluster_tags[mixture_n+i] = cluster.cluster_tag();
                sum = logadd(sum, (_lambda_log + _process_prior->log_predictive(cluster, _state) + log(_baseline_weights[i]) +
                                   cluster.model().log_predictive(range))/temp);
                log_weights[mixture_n+i] = sum;
        }
}

void
dpm_tfbs_t::mixture_weights(const vector<range_t>& range_set, double log_weights[], cluster_tag_t cluster_tags[], const double temp, const bool include_background)
{
        ssize_t mixture_n  = mixture_components();
        ssize_t baseline_n = baseline_components();
        double sum         = -HUGE_VAL;
        double n           = range_set.size();

        cluster_tag_t i = 0;
        ////////////////////////////////////////////////////////////////////////
        // loop through existing clusters
        for (cm_iterator it = _state.begin(); it != _state.end(); it++) {
                cluster_t& cluster = **it;
                cluster_tags[i] = cluster.cluster_tag();
                if (cluster.cluster_tag() == _state.bg_cluster_tag) {
                        ////////////////////////////////////////////////////////
                        // mixture component 1: background model
                        if (include_background) {
                                sum = logadd(sum, (n*_lambda_inv_log + cluster.model().log_predictive(range_set))/temp);
                        }
                }
                else {
                        ////////////////////////////////////////////////////////
                        // mixture component 2: dirichlet process
                        sum = logadd(sum, (n*_lambda_log + n*_process_prior->log_predictive(cluster, _state) + cluster.model().log_predictive(range_set))/temp);
                }
                log_weights[i] = sum;
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (i = 0; i < baseline_n; i++) {
                cluster_t& cluster = _state.get_free_cluster(_baseline_tags[i]);
                cluster_tags[mixture_n+i] = cluster.cluster_tag();
                sum = logadd(sum, (n*_lambda_log + n*_process_prior->log_predictive(cluster, _state) + log(_baseline_weights[i]) +
                                   cluster.model().log_predictive(range_set))/temp);
                log_weights[mixture_n+i] = sum;
        }
}

double
dpm_tfbs_t::likelihood() const {
        double result = 0;

        for (cm_iterator it = _state.begin();
             it != _state.end(); it++) {
                const cluster_t& cluster = **it;
                result += cluster.model().log_likelihood();
        }
        assert(!isnan(result));

        return result;
}

/*
 * compute P(X | Z) P(Z) \propto P(Z | X)
 */
double
dpm_tfbs_t::posterior() const {
        double result = likelihood();

        // background prior
        result += _state[_state.bg_cluster_tag].size()*_lambda_inv_log;
        // tfbs prior
        result += _state.num_tfbs*_lambda_log;
        // process prior
        result += _process_prior->joint(_state);

        assert(!isnan(result));

        return result;
}

void
dpm_tfbs_t::update_map()
{
        double posterior_value = posterior();

        if (_samples.map_value < posterior_value) {
                _samples.map_partition = _state.dpm_partition();
                _samples.map_value     = posterior_value;
        }
}

void
dpm_tfbs_t::record_partition()
{
        _samples.partitions.push_back(_state.dpm_partition());
}

void
dpm_tfbs_t::update_samples(size_t sampling_steps)
{
        // record for every position the average number of times it
        // belonged to the background model
        for (da_iterator it = _data.begin();
             it != _data.end(); it++) {
                const index_i& index  = **it;
                const size_t sequence = index[0];
                const size_t position = index[1];
                if (_state[index] == _state.bg_cluster_tag) {
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
        // update map partition and value
        update_map();
        // record current partition
        record_partition();
}

samples_t&
dpm_tfbs_t::samples() {
        return _samples;
}

dpm_tfbs_state_t&
dpm_tfbs_t::state() {
        return _state;
}

const dpm_tfbs_state_t&
dpm_tfbs_t::state() const {
        return _state;
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const dpm_tfbs_t& dpm)
{
        o << "Cluster Assignments:"           << endl;
        o << dpm._state.cluster_assignments() << endl;

        o << "TFBS Start Positions:"          << endl;
        o << dpm._state.tfbs_start_positions  << endl;

        o << "Clusters:"                      << endl;
        o << dpm._state;

        return o;
}
