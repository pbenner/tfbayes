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

#include <dpm-tfbs.hh>
#include <statistics.hh>

#include <tfbayes/logarithmetic.h>
#include <tfbayes/fastlog.h>

using namespace std;

DPM_TFBS::DPM_TFBS(
        double alpha, double d, double lambda, size_t tfbs_length,
        const DataTFBS& data, const DataTFBS& data_comp,
        gsl_matrix *baseline_priors[])
        : // length of tfbs
          TFBS_LENGTH(tfbs_length),
          // raw sequences
          _data(data),
          _data_comp(data_comp),
          // cluster manager
          _cluster_assignments(_data.lengths(), -1),
          _clustermanager(_cluster_assignments),
          // strength parameter for the dirichlet process
          alpha(alpha),
          alpha_log(log(alpha)),
          // pitman-yor discount factor
          d(d),
          d_log(log(d)),
          // mixture weight for the dirichlet process
          lambda(lambda),
          lambda_log(log(lambda)),
          lambda_inv_log(log(1-lambda)),
          // starting positions of tfbs
          _tfbs_start_positions(_data.lengths(), 0),
          // number of transcription factor binding sites
          num_tfbs(0)
{
        gsl_matrix* bg_alpha = init_alpha(BG_LENGTH);

        ////////////////////////////////////////////////////////////////////////////////
        // initialize joint posterior
        for (size_t i = 0; i < data.length(); i++) {
                _posterior.probabilities.push_back(vector<double>(data.length(i), 0.0));
        }

        ////////////////////////////////////////////////////////////////////////////////
        // add background model to the clustermanager
        bg_cluster_tag = _clustermanager.add_cluster(new ProductDirichlet(bg_alpha, _data));
        // add model components for the baseline measure
        for (size_t i = 0; baseline_priors[i] != NULL; i++) {
                const model_tag_t model_tag = _clustermanager.add_baseline_model(new ProductDirichlet(baseline_priors[i], _data));
                _model_tags.push_back(model_tag);
        }

        ////////////////////////////////////////////////////////////////////////////////
        // assign all elements to the background
        for (DataTFBS::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                _clustermanager[bg_cluster_tag].add_observations(range_t(*it, *it));
        }

        ////////////////////////////////////////////////////////////////////////////////
        // free prior
        gsl_matrix_free(bg_alpha);
}

DPM_TFBS::~DPM_TFBS() {
}

DPM_TFBS*
DPM_TFBS::clone() const {
        return new DPM_TFBS(*this);
}

bool
DPM_TFBS::valid_for_sampling(const index_t& index) const
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
                        if (_tfbs_start_positions[index_t(sequence, position+i)] == 1) {
                                return false;
                        }
                }
        }

        return true;
}

void
DPM_TFBS::add(const index_t& index, cluster_tag_t tag)
{
        const index_t& from(index);
        const index_t  to  (index[0], index[1] + TFBS_LENGTH - 1);

        if (tag == bg_cluster_tag) {
                _clustermanager[tag].add_observations(range_t(from, to));
        }
        else {
                _clustermanager[tag].add_observations(range_t(from, to));
                num_tfbs++;
                _tfbs_start_positions[index] = 1;
        }
}

void
DPM_TFBS::remove(const index_t& index, cluster_tag_t tag)
{
        const index_t& from(index);
        const index_t  to  (index[0], index[1] + TFBS_LENGTH - 1);

        if (tag == bg_cluster_tag) {
                _clustermanager[tag].remove_observations(range_t(from, to));
        }
        else {
                _clustermanager[tag].remove_observations(range_t(from, to));
                num_tfbs--;
                _tfbs_start_positions[index] = 0;
        }
}

size_t
DPM_TFBS::mixture_components() const
{
        return _clustermanager.size();
}

size_t
DPM_TFBS::baseline_components() const
{
        return _model_tags.size();
}

void
DPM_TFBS::mixture_weights(const index_t& index, double log_weights[], cluster_tag_t cluster_tags[])
{
        range_t range(index, index_t(index[0], index[1] + TFBS_LENGTH - 1));
        size_t mixture_n   = mixture_components();
        size_t baseline_n  = baseline_components();
        double dp_norm_log = log(num_tfbs + alpha);
//        double dp_norm_log = log(alpha + mixture_n - 1);
        double sum         = -HUGE_VAL;

        cluster_tag_t i = 0;
        for (ClusterManager::const_iterator it = _clustermanager.begin(); it != _clustermanager.end(); it++) {
                Cluster& cluster = **it;
                cluster_tags[i] = cluster.cluster_tag();
                ////////////////////////////////////////////////////////////////
                // mixture component 1: background model
                if (cluster_tags[i] == bg_cluster_tag) {
                        sum = logadd(sum, lambda_inv_log + cluster.model().log_pdf(range));
                        log_weights[i] = sum;
                }
                ////////////////////////////////////////////////////////////////
                // mixture component 2: dirichlet process for tfbs models
                else {
                        double num_elements = (double)cluster.size();
                        sum = logadd(sum, lambda_log + log(num_elements-d) - dp_norm_log + cluster.model().log_pdf(range));
//                        sum = logadd(sum, lambda_log - dp_norm_log + cluster.model().log_pdf(range));
                        log_weights[i] = sum;
                }
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        for (size_t i = 0; i < baseline_n; i++) {
                cluster_tags[mixture_n+i] = _clustermanager.get_free_cluster(_model_tags[i]).cluster_tag();
                sum = logadd(sum, lambda_log + log((alpha + d*(mixture_n-1))/baseline_n) - dp_norm_log +
                             _clustermanager[cluster_tags[mixture_n+i]].model().log_pdf(range));
//                sum = logadd(sum, lambda_log + log(alpha) - dp_norm_log + _clustermanager[cluster_tags[mixture_n+i]].model().log_pdf(range));
                log_weights[mixture_n+i] = sum;
        }
}

double
DPM_TFBS::likelihood() const {
        double result = 0;

        for (ClusterManager::const_iterator it = _clustermanager.begin();
             it != _clustermanager.end(); it++) {
                Cluster& cluster = **it;
                result += cluster.model().log_likelihood();
        }
        return result;
}

void
DPM_TFBS::update_graph(sequence_data_t<short> tfbs_start_positions)
{
        std::list<const index_t*> binding_sites;
        // find all binding sites
        for (Indexer::sampling_iterator it = _data.sampling_begin();
             it != _data.sampling_end(); it++) {
                if (tfbs_start_positions[**it] == 1) {
                        binding_sites.push_back(*it);
                }
        }
        // iterate over binding sites
        for (std::list<const index_t*>::const_iterator it = binding_sites.begin();
             it != binding_sites.end(); it++) {
                // if there still is a binding site
                if (tfbs_start_positions[**it] == 1) {
                        cluster_tag_t tag = _cluster_assignments[**it];
                        tfbs_start_positions[**it] = 0;
                        // find sites with the same cluster assignment
                        for (std::list<const index_t*>::const_iterator is = it;
                             is != binding_sites.end(); is++) {
                                if (tfbs_start_positions[**is] == 1 && _cluster_assignments[**is] == tag) {
                                        _posterior.graph[edge_t(**it, **is)]++;
                                }
                        }
                }
        }
}

void
DPM_TFBS::update_posterior(size_t sampling_steps) {
        if (sampling_steps % 100 == 0) {
                _tfbs_graph.cleanup(1);
        }
        for (DataTFBS::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                const index_t& index  = *it;
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

const posterior_t&
DPM_TFBS::posterior() const {
        return _posterior;
}

const DataTFBS&
DPM_TFBS::data() const {
        return _data;
}

const ClusterManager&
DPM_TFBS::clustermanager() const {
        return _clustermanager;
}

const Graph&
DPM_TFBS::graph() const {
        return _tfbs_graph;
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const DPM_TFBS& dpm)
{
        o << "Cluster Assignments:"     << endl;
        o << dpm._cluster_assignments   << endl;

        o << "TFBS Start Positions:"    << endl;
        o << dpm._tfbs_start_positions  << endl;

        o << "Clusters:"                << endl;
        o << dpm._clustermanager;

        return o;
}
