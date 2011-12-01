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
#include <tfbayes/fastlog.h>

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
          _cluster_assignments(_data.sizes(), -1),
          _clustermanager(_cluster_assignments),
          // strength parameter for the dirichlet process
          alpha(options.alpha),
          alpha_log(log(options.alpha)),
          // pitman-yor discount factor
          discount(options.discount),
          discount_log(log(options.discount)),
          // mixture weight for the dirichlet process
          lambda(options.lambda),
          lambda_log(log(options.lambda)),
          lambda_inv_log(log(1-options.lambda)),
          // starting positions of tfbs
          _tfbs_start_positions(_data.sizes(), 0),
          // number of transcription factor binding sites
          num_tfbs(0)
{
        ////////////////////////////////////////////////////////////////////////////////
        // initialize joint posterior
        for (size_t i = 0; i < data.size(); i++) {
                _posterior.probabilities.push_back(vector<double>(data.size(i), 0.0));
        }

        ////////////////////////////////////////////////////////////////////////////////
        // add background model to the clustermanager
        if (options.background_model == "independence" || options.background_model == "") {
                const matrix<double>& bg_alpha = init_alpha(BG_LENGTH);
                ProductDirichlet* bg = new ProductDirichlet(bg_alpha, _data);
                bg_cluster_tag = _clustermanager.add_cluster(bg);
        }
        else if (options.background_model == "markov chain mixture") {
                assert(options.context >= 0);
                MarkovChainMixture* bg = new MarkovChainMixture(ALPHABET_SIZE, options.context, _data, _cluster_assignments, 0);
                bg_cluster_tag = _clustermanager.add_cluster(bg);
        }
        else if (options.background_model == "parsimonious tree") {
                assert(options.context >= 0);
                ParsimoniousTree* bg = new ParsimoniousTree(ALPHABET_SIZE, options.context, _data, _cluster_assignments, 0);
                bg_cluster_tag = _clustermanager.add_cluster(bg);
        }
        else {
                cerr << "Unknown background model." << endl;
                exit(EXIT_FAILURE);
        }
        // add model components for the baseline measure
        for (size_t i = 0; i < options.baseline_priors.size(); i++) {
                assert(options.tfbs_length == options.baseline_priors[i].size());
                ProductDirichlet* dirichlet = new ProductDirichlet(options.baseline_priors[i], _data);
                const model_tag_t model_tag = _clustermanager.add_baseline_model(dirichlet);
                _model_tags.push_back(model_tag);
        }

        ////////////////////////////////////////////////////////////////////////////////
        // assign all elements to the background
        for (data_tfbs_t::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                range_t range(**it, 1);
                _clustermanager[bg_cluster_tag].add_observations(range);
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
}

DpmTfbs::~DpmTfbs() {
}

void
DpmTfbs::test_background() {
        seq_index_t index1(0,0);
        seq_index_t index2(0,10);
        seq_index_t index3(0,11);
        seq_index_t index4(0,12);
        seq_index_t index5(0,13);
        seq_index_t index6(0,22);

        Cluster& cluster1 = _clustermanager.get_free_cluster(_model_tags[0]);
        cluster_tag_t cluster_tag1 = cluster1.cluster_tag();
        cout << "Adding index1:" << index1 << " to cluster:" << cluster_tag1 << endl;
        remove(index1, bg_cluster_tag);
        add(index1, cluster_tag1);
        cout << "Adding index6:" << index6<< " to cluster:" << cluster_tag1 << endl;
        remove(index6, bg_cluster_tag);
        add(index6, cluster_tag1);
        cout << endl;

        cout << _cluster_assignments << endl;

        Cluster& cluster2 = _clustermanager.get_free_cluster(_model_tags[0]);
        cluster_tag_t cluster_tag2 = cluster2.cluster_tag();

        cout << "Adding index2:" << index2 << " to cluster:" << cluster_tag2 << endl;
        remove(index2, bg_cluster_tag);
        add(index2, cluster_tag2);
        cout << _cluster_assignments;
        remove(index2, cluster_tag2);
        add(index2, bg_cluster_tag);
        cout << "Removing index2:" << index2 << " from cluster:" << cluster_tag2 << endl << endl;

        cout << "Adding index3:" << index3 << " to cluster:" << cluster_tag2 << endl;
        remove(index3, bg_cluster_tag);
        add(index3, cluster_tag2);
        cout << _cluster_assignments;
        remove(index3, cluster_tag2);
        add(index3, bg_cluster_tag);
        cout << "Removing index3:" << index2 << " from cluster:" << cluster_tag2 << endl << endl;

        cout << "Adding index4:" << index4 << " to cluster:" << cluster_tag2 << endl;
        remove(index4, bg_cluster_tag);
        add(index4, cluster_tag2);
        cout << _cluster_assignments;
        remove(index4, cluster_tag2);
        add(index4, bg_cluster_tag);
        cout << "Removing index4:" << index2 << " from cluster:" << cluster_tag2 << endl << endl;

        cout << "Adding index5:" << index4 << " to cluster:" << cluster_tag2 << endl;
        remove(index5, bg_cluster_tag);
        add(index5, cluster_tag2);
        cout << _cluster_assignments;
        remove(index5, cluster_tag2);
        add(index5, bg_cluster_tag);
        cout << "Removing index5:" << index2 << " from cluster:" << cluster_tag2 << endl;

        exit(EXIT_SUCCESS);
}

void
DpmTfbs::test() {
        seq_index_t index1(0,0);
        seq_index_t index2(1,0);
        seq_index_t index3(2,0);
        seq_index_t index4(3,0);

        Cluster& cluster1 = _clustermanager.get_free_cluster(_model_tags[0]);
        cluster_tag_t cluster_tag1 = cluster1.cluster_tag();
        cout << "Adding index1:" << index1 << " to cluster:" << cluster_tag1 << endl;
        add(index1, cluster_tag1);
        Cluster& cluster2 = _clustermanager.get_free_cluster(_model_tags[0]);
        cluster_tag_t cluster_tag2 = cluster2.cluster_tag();
        cout << "Adding index2:" << index2 << " to cluster:" << cluster_tag2 << endl;
        add(index2, cluster_tag2);

        cout << "Components: " << mixture_components() << " + " << baseline_components() << endl;
        size_t components = mixture_components() + baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        cluster_tag_t new_cluster_tag;

        cout << "Sampling index3:" << index3 << endl;
        mixture_weights(index3, log_weights, cluster_tags);
        for (size_t i = 0; i < 100; i++) {
                new_cluster_tag = cluster_tags[select_component(components, log_weights)];
                cout << "selected cluster " << new_cluster_tag << endl;
        }

        cout << "Sampling index4:" << index4 << endl;
        mixture_weights(index4, log_weights, cluster_tags);
        for (size_t i = 0; i < 100; i++) {
                new_cluster_tag = cluster_tags[select_component(components, log_weights)];
                cout << "selected cluster " << new_cluster_tag << endl;
        }

        exit(EXIT_SUCCESS);
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
        if (_tfbs_start_positions[index] == 0) {
                // check if this element belongs to a tfbs that starts
                // earlier in the sequence
                if (_clustermanager[index] != bg_cluster_tag) {
                        return false;
                }
                // check if there is a tfbs starting within the word
                for (size_t i = 1; i < TFBS_LENGTH; i++) {
                        if (_tfbs_start_positions[seq_index_t(sequence, position+i)] == 1) {
                                return false;
                        }
                }
        }

        return true;
}

void
DpmTfbs::add(const index_i& index, cluster_tag_t tag)
{
        const range_t range(index, TFBS_LENGTH);

        if (tag == bg_cluster_tag) {
                _clustermanager[tag].add_observations(range);
        }
        else {
                _clustermanager[tag].add_observations(range);
                num_tfbs++;
                _tfbs_start_positions[index] = 1;
        }
}

void
DpmTfbs::remove(const index_i& index, cluster_tag_t tag)
{
        const range_t range(index, TFBS_LENGTH);

        if (tag == bg_cluster_tag) {
                _clustermanager[tag].remove_observations(range);
        }
        else {
                _clustermanager[tag].remove_observations(range);
                num_tfbs--;
                _tfbs_start_positions[index] = 0;
        }
}

size_t
DpmTfbs::mixture_components() const
{
        return _clustermanager.size();
}

size_t
DpmTfbs::baseline_components() const
{
        return _model_tags.size();
}

double
DpmTfbs::py_prior(Cluster& cluster)
{
        if (cluster.size() == 0) {
                return log(alpha + discount*(mixture_components()-1)) - log(num_tfbs + alpha);
        }
        else {
                return log(cluster.size()-discount) - log(num_tfbs + alpha);
        }
}

double
DpmTfbs::uniform_prior(Cluster& cluster)
{
        double K = mixture_components()-1;

        if (cluster.size() == 0) {
                return log(alpha) - log(alpha + K);
        }
        else {
                return -log(alpha + K);
        }
}

double
DpmTfbs::poppe_prior(Cluster& cluster)
{
        double K = mixture_components()-1;
        double N = num_tfbs;

        if (K == 0 && cluster.size() == 0) {
                return 0;
        }
        if (cluster.size() == 0) {
                if (K == 1.0) {
                        return -log(N+1.0);
                }
                else {
                        return log(K*(K-1.0)/(N*(N+1.0)));
                }
        }
        else {
                if (K == 1.0) {
                        return log(N/(N+1.0));
                }
                else {
                        return log((cluster.size()+1.0)/(N+1.0) * (N-K+1.0)/N);
                }
        }
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
        for (ClusterManager::const_iterator it = _clustermanager.begin(); it != _clustermanager.end(); it++) {
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
                Cluster& cluster = _clustermanager.get_free_cluster(_model_tags[i]);
                cluster_tags[mixture_n+i] = cluster.cluster_tag();
                sum = logadd(sum, lambda_log + __process_prior(cluster) + log(_baseline_weights[i]) +
                             cluster.model().log_predictive(range));
                log_weights[mixture_n+i] = sum;
        }
}

double
DpmTfbs::likelihood() const {
        double result = 0;

        for (ClusterManager::const_iterator it = _clustermanager.begin();
             it != _clustermanager.end(); it++) {
                Cluster& cluster = **it;
                result += cluster.model().log_likelihood();
        }
        return result;
}

void
DpmTfbs::update_graph(sequence_data_t<short> tfbs_start_positions)
{
        std::list<const index_i*> binding_sites;
        // find all binding sites
        for (Indexer::sampling_iterator it = _data.sampling_begin();
             it != _data.sampling_end(); it++) {
                if (tfbs_start_positions[**it] == 1) {
                        binding_sites.push_back(*it);
                }
        }
        // iterate over binding sites
        for (std::list<const index_i*>::const_iterator it = binding_sites.begin();
             it != binding_sites.end(); it++) {
                // if there still is a binding site
                if (tfbs_start_positions[**it] == 1) {
                        cluster_tag_t tag = _cluster_assignments[**it];
                        tfbs_start_positions[**it] = 0;
                        // find sites with the same cluster assignment
                        for (std::list<const index_i*>::const_iterator is = it;
                             is != binding_sites.end(); is++) {
                                if (tfbs_start_positions[**is] == 1 && _cluster_assignments[**is] == tag) {
                                        _posterior.graph[edge_t(**it, **is)]++;
                                }
                        }
                }
        }
}

void
DpmTfbs::update_hypergraph(sequence_data_t<short> tfbs_start_positions)
{
        list<const index_i*> binding_sites;
        stringstream ss;
        // find all binding sites
        for (Indexer::sampling_iterator it = _data.sampling_begin();
             it != _data.sampling_end(); it++) {
                if (tfbs_start_positions[**it] == 1) {
                        binding_sites.push_back(*it);
                }
        }
        // iterate over binding sites
        for (list<const index_i*>::const_iterator it = binding_sites.begin();
             it != binding_sites.end(); it++) {
                // if there still is a binding site
                if (tfbs_start_positions[**it] == 1) {
                        cluster_tag_t tag = _cluster_assignments[**it];
                        tfbs_start_positions[**it] = 0;
                        ss << "{ " << *static_cast<const seq_index_t*>(*it);
                        // find sites with the same cluster assignment
                        for (list<const index_i*>::const_iterator is = it;
                             is != binding_sites.end(); is++) {
                                if (tfbs_start_positions[**is] == 1 && _cluster_assignments[**is] == tag) {
                                        ss << " " << *static_cast<const seq_index_t*>(*is);
                                }
                        }
                        ss << " } ";
                }
        }
        if (!ss.str().size() == 0) {
                _posterior.hypergraph.push_back(ss.str());
        }
}

void
DpmTfbs::update_posterior(size_t sampling_steps) {
        if (sampling_steps % 100 == 0) {
                _tfbs_graph.cleanup(1);
        }
        for (data_tfbs_t::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                const index_i& index  = **it;
                const size_t sequence = index[0];
                const size_t position = index[1];
                if (_clustermanager[index] == bg_cluster_tag) {
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
        update_graph(_tfbs_start_positions);
}

posterior_t&
DpmTfbs::posterior() {
        return _posterior;
}

const data_tfbs_t&
DpmTfbs::data() const {
        return _data;
}

const ClusterManager&
DpmTfbs::clustermanager() const {
        return _clustermanager;
}

const Graph&
DpmTfbs::graph() const {
        return _tfbs_graph;
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
        o << dpm._cluster_assignments   << endl;

        o << "TFBS Start Positions:"    << endl;
        o << dpm._tfbs_start_positions  << endl;

        o << "Clusters:"                << endl;
        o << dpm._clustermanager;

        return o;
}
