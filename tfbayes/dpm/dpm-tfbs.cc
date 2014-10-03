/* Copyright (C) 2011-2013 Philipp Benner
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
#include <limits>

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/utility/statistics.hh>
#include <tfbayes/utility/logarithmetic.hh>

using namespace std;

dpm_tfbs_t::dpm_tfbs_t(const tfbs_options_t& options,
                       const data_tfbs_t& data,
                       boost::optional<const alignment_set_t<>&> alignment_set)
        : // baseline
          _baseline_weights(options.baseline_weights),
          // phylogenetic information
          _data(&data),
          // coded nucleotide sequences
          _alignment_set(NULL),
          // cluster manager
          _state(data.sizes(), options.tfbs_length, 0, data),
          // mixture weight for the dirichlet process
          _lambda(options.lambda),
          _lambda_log(log(options.lambda)),
          _lambda_inv_log(log(1-options.lambda)),
          // length of tfbs
          _tfbs_length(options.tfbs_length)
{
        ////////////////////////////////////////////////////////////////////////////////
        // check that the alignment data matches the phylogenetic data
        if (alignment_set) {
                assert(data.size() == alignment_set->size());
                for (size_t i = 0; i < data.size(); i++) {
                        assert(data[i].size() == (*alignment_set)[i].length());
                }
                _alignment_set = &*alignment_set;
        }

        ////////////////////////////////////////////////////////////////////////////////
        // add background model to the state
        if (options.background_model == "independence-dirichlet" || options.background_model == "") {
                /* every position in the background is fully
                 * independet, this give more flexibility to the
                 * prior pseudocounts */
                independence_background_t* bg = new independence_background_t(options.background_alpha, data, _state.cluster_assignments());
                _state.bg_cluster_tag = _state.add_cluster(bg);
                bg->set_bg_cluster_tag(_state.bg_cluster_tag);
        }
        else if (options.background_model == "independence-dirichlet-gamma") {
                /* every position in the background is fully
                 * independet, where the Dirichlet pseudocounts
                 * are integrated out */
                cerr << "TODO: variable gamma prior parameters" << endl;
                independence_background_t* bg = new independence_background_t(2.0, 2.0, data, _state.cluster_assignments());
                _state.bg_cluster_tag = _state.add_cluster(bg);
                bg->set_bg_cluster_tag(_state.bg_cluster_tag);
        }
        else if (options.background_model == "dirichlet") {
                /* single dirichlet-compound distribution for all
                 * nucleotides in the background */
                product_dirichlet_t* bg = new product_dirichlet_t(options.background_alpha, data, data.complements());
                _state.bg_cluster_tag = _state.add_cluster(bg);
        }
        else if (options.background_model == "markov chain mixture") {
                assert(options.background_context >= 0);
                markov_chain_mixture_t* bg = new markov_chain_mixture_t(data_tfbs_t::alphabet_size, options, data, _state.cluster_assignments(), 0);
                _state.bg_cluster_tag = _state.add_cluster(bg);
        }
        else {
                cerr << "Unknown background model." << endl;
                exit(EXIT_FAILURE);
        }
        // baseline weights are already initialized
        baseline_priors_t::const_iterator it = options.baseline_priors.begin();
        baseline_tags_t  ::const_iterator is = options.baseline_tags  .begin();
        for (;
             it != options.baseline_priors.end() &&
             is != options.baseline_tags  .end(); it++, is++) {
                assert(it->size() == options.tfbs_length);
                for (size_t j = 0; j < options.tfbs_length; j++) {
                        assert((*it)[j].size() == data_tfbs_t::alphabet_size);
                }
                product_dirichlet_t* dirichlet = new product_dirichlet_t(*it, data, data.complements());
                _state.add_baseline_model(dirichlet, *is);
                _baseline_tags.push_back(*is);
        }
        ////////////////////////////////////////////////////////////////////////////////
        // assign all elements to the background
        for (da_iterator it = data.begin();
             it != data.end(); it++) {
                // cannot use _state.add() here because we need to add
                // sites of length one to the background
                range_t range(**it, 1);
                _state[_state.bg_cluster_tag].add_observations(range);
        }
        ////////////////////////////////////////////////////////////////////////////////
        // set the process prior
        if (options.process_prior == "pitman-yor process" || options.process_prior == "") {
                _process_prior = new pitman_yor_prior(options.alpha, options.discount, _state.bg_cluster_tag);
        }
        else if (options.process_prior == "uniform process") {
                _process_prior = new uniform_prior(options.alpha);
        }
        else if (options.process_prior == "poppe process") {
                _process_prior = new poppe_prior();
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

void
swap(dpm_tfbs_t& first, dpm_tfbs_t& second)
{
        swap(first._baseline_weights, second._baseline_weights);
        swap(first._baseline_tags,    second._baseline_tags);
        swap(first._data,             second._data);
        swap(first._alignment_set,    second._alignment_set);
        swap(first._state,            second._state);
        swap(first._lambda,           second._lambda);
        swap(first._lambda_log,       second._lambda_log);
        swap(first._lambda_inv_log,   second._lambda_inv_log);
        swap(first._tfbs_length,      second._tfbs_length);
        swap(first._process_prior,    second._process_prior);
}

dpm_tfbs_t&
dpm_tfbs_t::operator=(const mixture_model_t& mixture_model)
{
        dpm_tfbs_t tmp(static_cast<const dpm_tfbs_t&>(mixture_model));
        swap(*this, tmp);
        return *this;
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

GCC_ATTRIBUTE_HOT
void
dpm_tfbs_t::mixture_weights(const range_t& range, double log_weights[], cluster_tag_t cluster_tags[], const double temp, const double baseline)
{
        double sum = baseline;

        cluster_tag_t i = 0;
        ////////////////////////////////////////////////////////////////////////
        // loop through existing clusters
        for (cm_iterator it = _state.begin(); it != _state.end(); it++, i++) {
                cluster_t& cluster = **it;
                if (cluster.cluster_tag() == _state.bg_cluster_tag) {
                        ////////////////////////////////////////////////////////
                        // mixture component 1: background model
                        sum = logadd(sum, (_lambda_inv_log + cluster.model().log_predictive(range))/temp);
                        log_weights [i] = sum;
                        cluster_tags[i] = _state.bg_cluster_tag;
                }
                else {
                        ////////////////////////////////////////////////////////
                        // mixture component 2: dirichlet process
                        sum = logadd(sum, (_lambda_log + _process_prior->log_predictive(cluster, _state) + cluster.model().log_predictive(range))/temp);
                        log_weights [i] = sum;
                        cluster_tags[i] = cluster.cluster_tag();
                }
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (size_t j = 0; j < baseline_components(); j++, i++) {
                cluster_t& cluster = _state.get_free_cluster(_baseline_tags[j]);
                sum = logadd(sum, (_lambda_log + _process_prior->log_predictive(cluster, _state) + log(_baseline_weights[j]) +
                                   cluster.model().log_predictive(range))/temp);
                log_weights [i] = sum;
                cluster_tags[i] = cluster.cluster_tag();
        }
        assert(static_cast<size_t>(i) == mixture_components() + baseline_components());
}

GCC_ATTRIBUTE_HOT
void
dpm_tfbs_t::mixture_weights(const vector<range_t>& range_set, double log_weights[], cluster_tag_t cluster_tags[], const double temp, const bool include_background)
{
        double sum         = -numeric_limits<double>::infinity();
        double n           = range_set.size();

        cluster_tag_t i = 0;
        ////////////////////////////////////////////////////////////////////////
        // loop through existing clusters
        for (cm_iterator it = _state.begin(); it != _state.end(); it++, i++) {
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
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (size_t j = 0; j < baseline_components(); j++, i++) {
                cluster_t& cluster = _state.get_free_cluster(_baseline_tags[j]);
                sum = logadd(sum, (n*_lambda_log + n*_process_prior->log_predictive(cluster, _state) + log(_baseline_weights[j]) +
                                   cluster.model().log_predictive(range_set))/temp);
                log_weights [i] = sum;
                cluster_tags[i] = cluster.cluster_tag();
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
        assert(!std::isnan(result));

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

        assert(!std::isnan(result));

        return result;
}

dpm_tfbs_state_t&
dpm_tfbs_t::state() {
        return _state;
}

const dpm_tfbs_state_t&
dpm_tfbs_t::state() const {
        return _state;
}

const data_tfbs_t&
dpm_tfbs_t::data() const {
        return *_data;
}

const alignment_set_t<>&
dpm_tfbs_t::alignment_set() const {
        return *_alignment_set;
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
