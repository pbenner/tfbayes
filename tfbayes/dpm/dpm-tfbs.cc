/* Copyright (C) 2011-2014 Philipp Benner
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

#include <boost/foreach.hpp>

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/utility/statistics.hh>
#include <tfbayes/utility/logarithmetic.hh>
#include <tfbayes/utility/normalize.hh>

using namespace std;

dpm_tfbs_t::dpm_tfbs_t(const tfbs_options_t& options,
                       const data_tfbs_t& data,
                       boost::optional<const alignment_set_t<>&> alignment_set)
          // phylogenetic information
        : m_data           (&data)
          // coded nucleotide sequences
        , m_alignment_set  (NULL)
          // cluster manager
        , m_state          (options, data)
          // mixture weight for the dirichlet process
        , m_lambda         (options.lambda)
        , m_lambda_log     (log(options.lambda))
        , m_lambda_inv_log (log(1-options.lambda))
{
        ////////////////////////////////////////////////////////////////////////////////
        // check that the alignment data matches the phylogenetic data
        if (alignment_set) {
                assert(data.size() == alignment_set->size());
                for (size_t i = 0; i < data.size(); i++) {
                        assert(data[i].size() == (*alignment_set)[i].length());
                }
                m_alignment_set = &*alignment_set;
        }

        ////////////////////////////////////////////////////////////////////////////////
        // add background model to the state
        if (options.background_model == "independence-dirichlet") {
                /* every position in the background is fully
                 * independet; the Dirichlet pseudocounts
                 * can be integrated out */
                assert(options.background_gamma.size() == 2);
                assert(options.threads >= 1);
                assert(options.background_alpha.size() == 1);
                thread_pool_t thread_pool(options.threads);
                independence_background_t* bg = new independence_background_t(
                        options.background_alpha[0], options.background_gamma, data,
                        m_state.cluster_assignments(), thread_pool, options.background_cache,
                        alignment_set);
                cluster_tag_t tag = m_state.add_background_cluster(*bg);
                bg->set_bg_cluster_tag(tag);
        }
        else if (options.background_model == "independence-mixture-dirichlet") {
                independence_mixture_background_t* bg = new independence_mixture_background_t(
                        options.background_alpha, options.background_weights, data,
                        m_state.cluster_assignments(), alignment_set);
                cluster_tag_t tag = m_state.add_background_cluster(*bg);
                bg->set_bg_cluster_tag(tag);
        }
        else if (options.background_model == "default-background") {
                assert(options.background_gamma.size() == 2);
                assert(options.threads >= 1);
                thread_pool_t thread_pool(options.threads);
                default_background_t* bg = new default_background_t(
                        options.background_alpha,
                        options.background_gamma,
                        options.background_weights, data,
                        m_state.cluster_assignments(), thread_pool, options.background_cache,
                        alignment_set, options.verbose);
                cluster_tag_t tag = m_state.add_background_cluster(*bg);
                bg->set_bg_cluster_tag(tag);
        }
        else if (options.background_model == "entropy-background") {
                assert(options.background_beta.size() == 2);
                assert(options.threads >= 1);
                thread_pool_t thread_pool(options.threads);
                entropy_background_t* bg = new entropy_background_t(
                        options.background_beta, data,
                        m_state.cluster_assignments(), thread_pool, options.background_cache,
                        alignment_set, options.verbose);
                cluster_tag_t tag = m_state.add_background_cluster(*bg);
                bg->set_bg_cluster_tag(tag);
        }
        else if (options.background_model == "dirichlet-mixture") {
                /* multiple dirichlet-compound distribution for all
                 * nucleotides in the background */
                assert(options.background_alpha.size() >= 1);
                mixture_dirichlet_t* bg = new mixture_dirichlet_t(
                        options.background_alpha,
                        options.background_weights,
                        data);
                m_state.add_background_cluster(*bg);
        }
        else if (options.background_model == "markov chain mixture") {
                assert(options.background_context >= 0);
                markov_chain_mixture_t* bg = new markov_chain_mixture_t(
                        data_tfbs_t::alphabet_size, options, data, m_state.cluster_assignments(), 0);
                m_state.add_background_cluster(*bg);
        }
        else {
                cerr << "Unknown background model." << endl;
                exit(EXIT_FAILURE);
        }
        ////////////////////////////////////////////////////////////////////////////////
        // add the foreground models
        assert(options.baseline_priors.size() > 0);
        assert(options.baseline_priors.size() == options.baseline_lengths.size());
        assert(options.baseline_priors.size() == options.baseline_names  .size());
        assert(options.baseline_priors.size() == options.baseline_weights.size());
        baseline_priors_t ::const_iterator it = options.baseline_priors .begin();
        baseline_names_t  ::const_iterator is = options.baseline_names  .begin();
        baseline_weights_t::const_iterator ir = options.baseline_weights.begin();
        baseline_lengths_t::const_iterator iq = options.baseline_lengths.begin();
        for (; it != options.baseline_priors.end(); it++, is++, ir++, iq++) {
                assert(it->size() == 1);
                // add a baseline_prior for each foreground length
                BOOST_FOREACH(const double& length, *iq) {
                        if (options.verbose >= 1) {
                                flockfile(stderr);
                                cerr << boost::format("Adding baseline model `%s' of length %d.")
                                        % *is % length << endl;
                                fflush(stderr);
                                funlockfile(stderr);
                        }
                        model_id_t model_id = {*is, size_t(length)};
                        product_dirichlet_t* dirichlet = new product_dirichlet_t(model_id, *it, *iq, data, data.complements());
                        m_baseline_tags   .push_back(m_state.add_baseline_model(dirichlet));
                        // the weight consists of the individual
                        // weight of the baseline component and a
                        // normalization for the number of lenths
                        m_baseline_weights.push_back(log(*ir) + log(1.0/iq->size()));
                        // assume that the baseline tag is simply the
                        // size of the vector, so that we can easily
                        // index the baseline weights
                        assert(m_baseline_tags.back() == baseline_tag_t(m_baseline_tags.size())-1);
                }
        }
        // the baseline weights have to be normalized!
        m_baseline_weights = log_normalize(m_baseline_weights);
        ////////////////////////////////////////////////////////////////////////////////
        // assign all elements to the background
        for (da_iterator it = data.begin();
             it != data.end(); it++) {
                // cannot use m_state.add() here because we need to add
                // sites of length one to the background
                range_t range(*it, 1);
                m_state[*m_state.bg_cluster_tags.begin()].add_observations(range);
        }
        m_state[*m_state.bg_cluster_tags.begin()].update();
        ////////////////////////////////////////////////////////////////////////////////
        // set the process prior
        if (options.process_prior == "pitman-yor process" || options.process_prior == "") {
                m_process_prior = new pitman_yor_prior(options.alpha, options.discount);
        }
        else if (options.process_prior == "uniform process") {
                m_process_prior = new uniform_prior(options.alpha);
        }
        else if (options.process_prior == "poppe process") {
                m_process_prior = new poppe_prior();
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
        : m_baseline_weights (dpm.m_baseline_weights)
        , m_baseline_tags    (dpm.m_baseline_tags)
        // phylogenetic data
        , m_data             (dpm.m_data)
        // coded nucleotide sequences
        , m_alignment_set    (dpm.m_alignment_set)
        // cluster manager
        , m_state            (dpm.m_state)
        // mixture weight for the dirichlet process
        , m_lambda           (dpm.m_lambda)
        , m_lambda_log       (dpm.m_lambda_log)
        , m_lambda_inv_log   (dpm.m_lambda_inv_log)
        // process prios
        , m_process_prior    (dpm.m_process_prior->clone())
{ }

dpm_tfbs_t::~dpm_tfbs_t() {
        delete(m_process_prior);
}

dpm_tfbs_t*
dpm_tfbs_t::clone() const {
        return new dpm_tfbs_t(*this);
}

void
swap(dpm_tfbs_t& first, dpm_tfbs_t& second)
{
        swap(first.m_baseline_weights, second.m_baseline_weights);
        swap(first.m_baseline_tags,    second.m_baseline_tags);
        swap(first.m_data,             second.m_data);
        swap(first.m_alignment_set,    second.m_alignment_set);
        swap(first.m_state,            second.m_state);
        swap(first.m_lambda,           second.m_lambda);
        swap(first.m_lambda_log,       second.m_lambda_log);
        swap(first.m_lambda_inv_log,   second.m_lambda_inv_log);
        swap(first.m_process_prior,    second.m_process_prior);
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
        return m_state.size();
}

size_t
dpm_tfbs_t::baseline_components() const
{
        return m_baseline_tags.size();
}

double
dpm_tfbs_t::background_mixture_weight(const range_t& range, cluster_t& cluster)
{
        return range.length()*m_lambda_inv_log + cluster.model().log_predictive(range);
}

double
dpm_tfbs_t::foreground_mixture_weight(const range_t& range, cluster_t& cluster)
{
        double result = m_lambda_log + m_baseline_weights[cluster.baseline_tag()] + m_process_prior->log_predictive(cluster, m_state);
        // get the length of the foreground model
        size_t cluster_length = cluster.model().id().length;
        if (cluster_length < range.length()) {
                cluster_t& bg_cluster = m_state[*m_state.bg_cluster_tags.begin()];
                // split range in two
                range_t range1(range);
                range_t range2(range);
                range1.length()     = cluster_length;
                range2.index ()[1] += cluster_length;
                range2.length()    -= cluster_length;
                // compute log_predictives
                result += cluster.model().log_predictive(range1);
                // remaining positions are assigned to the background
                // model
                result += background_mixture_weight(range2, bg_cluster);
        }
        else {
                result += cluster.model().log_predictive(range);
        }
        return result;
}

GCC_ATTRIBUTE_HOT
void
dpm_tfbs_t::mixture_weights(
        const range_t& range, double log_weights[], cluster_tag_t cluster_tags[],
        const double temp, const bool include_background, const double baseline)
{
        double sum = baseline;

        cluster_tag_t i = 0;
        ////////////////////////////////////////////////////////////////////////
        // loop through existing clusters
        for (cm_iterator it = m_state.begin(); it != m_state.end(); it++, i++) {
                cluster_t& cluster = **it;
#ifdef DEBUG
                test_posterior(cluster, range);
#endif
                if (m_state.is_background(cluster)) {
                        ////////////////////////////////////////////////////////
                        // mixture component 1: background model
                        if (include_background) {
                                sum = logadd(sum, background_mixture_weight(range, cluster)/temp);
                        }
                }
                else {
                        ////////////////////////////////////////////////////////
                        // mixture component 2: dirichlet process
                        sum = logadd(sum, foreground_mixture_weight(range, cluster)/temp);
                }
                log_weights [i] = sum;
                cluster_tags[i] = cluster.cluster_tag();
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (size_t j = 0; j < baseline_components(); j++, i++) {
                cluster_t& cluster = m_state.get_free_cluster(m_baseline_tags[j]);
                sum = logadd(sum, foreground_mixture_weight(range, cluster)/temp);
#ifdef DEBUG
                test_posterior(cluster, range);
#endif
                log_weights [i] = sum;
                cluster_tags[i] = cluster.cluster_tag();
        }
        assert(static_cast<size_t>(i) == mixture_components() + baseline_components());
}

double
dpm_tfbs_t::background_mixture_weight(const vector<range_t>& range_set, cluster_t& cluster)
{
        double result = 0.0;

        // background weight
        result += range_set.size()*range_set[0].length()*m_lambda_inv_log;
        // predictive distribution
        result += cluster.model().log_predictive(range_set);

        return result;
}

double
dpm_tfbs_t::foreground_mixture_weight(const vector<range_t>& range_set, cluster_t& cluster)
{
        double result = 0.0;

        // foreground weight
        result += range_set.size()*m_lambda_log;
        // baseline weight
        result += range_set.size()*m_baseline_weights[cluster.baseline_tag()];
        // process prior
        result += m_process_prior->log_predictive(cluster, m_state, range_set.size());
        // predictive distribution
        result += cluster.model().log_predictive(range_set);

        return result;
}

GCC_ATTRIBUTE_HOT
void
dpm_tfbs_t::mixture_weights(const vector<range_t>& range_set, double log_weights[], cluster_tag_t cluster_tags[],
                            const double temp, const bool include_background)
{
        double sum = -numeric_limits<double>::infinity();

        cluster_tag_t i = 0;
        ////////////////////////////////////////////////////////////////////////
        // loop through existing clusters
        for (cm_iterator it = m_state.begin(); it != m_state.end(); it++, i++) {
                cluster_t& cluster = **it;
#ifdef DEBUG
                test_posterior(cluster, range_set);
#endif
                if (m_state.is_background(cluster)) {
                        ////////////////////////////////////////////////////////
                        // mixture component 1: background model
                        if (include_background) {
                                sum = logadd(sum, background_mixture_weight(range_set, cluster)/temp);
                        }
                }
                else {
                        ////////////////////////////////////////////////////////
                        // mixture component 2: dirichlet process
                        sum = logadd(sum, foreground_mixture_weight(range_set, cluster)/temp);
                }
                log_weights [i] = sum;
                cluster_tags[i] = cluster.cluster_tag();
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (size_t j = 0; j < baseline_components(); j++, i++) {
                cluster_t& cluster = m_state.get_free_cluster(m_baseline_tags[j]);
#ifdef DEBUG
                test_posterior(cluster, range_set);
#endif
                sum = logadd(sum, foreground_mixture_weight(range_set, cluster)/temp);
                log_weights [i] = sum;
                cluster_tags[i] = cluster.cluster_tag();
        }
}

double
dpm_tfbs_t::likelihood() const {
        double result = 0;

        BOOST_FOREACH (const cluster_t* cluster, m_state) {
                result += cluster->model().log_likelihood();
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

        BOOST_FOREACH (const cluster_t* cluster, m_state) {
                if (m_state.is_background(*cluster)) {
                        // background prior
                        result += cluster->size()*m_lambda_inv_log;
                }
                else {
                        // foreground weight
                        result += cluster->size()*m_lambda_log;
                        result += cluster->size()*m_baseline_weights[cluster->baseline_tag()];
                }
        }
        // process prior
        result += m_process_prior->joint(m_state);

        assert(!std::isnan(result));

        return result;
}

dpm_tfbs_state_t&
dpm_tfbs_t::state() {
        return m_state;
}

const dpm_tfbs_state_t&
dpm_tfbs_t::state() const {
        return m_state;
}

const data_tfbs_t&
dpm_tfbs_t::data() const {
        return *m_data;
}

const alignment_set_t<>&
dpm_tfbs_t::alignment_set() const {
        return *m_alignment_set;
}

const dpm_tfbs_prior_t&
dpm_tfbs_t::process_prior() const {
        return *m_process_prior;
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const dpm_tfbs_t& dpm)
{
        o << "Cluster Assignments:"           << endl;
        o << dpm.m_state.cluster_assignments() << endl;

        o << "TFBS Start Positions:"          << endl;
        o << dpm.m_state.tfbs_start_positions  << endl;

        o << "Clusters:"                      << endl;
        o << dpm.m_state;

        return o;
}
