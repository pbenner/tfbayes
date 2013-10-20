/* Copyright (C) 2013 Philipp Benner
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

#define __STDC_LIMIT_MACROS

#include <stdint.h>

#include <boost/unordered/unordered_set.hpp>

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/dpm/dpm-sampling-history.hh>
#include <tfbayes/utility/progress.hh>

using namespace std;

// utility
// -----------------------------------------------------------------------------

size_t hash_value(const seq_index_t& seq_index)
{
        return seq_index.hash();
}

// functions for computing means and medians
// -----------------------------------------------------------------------------

static void
init_data(const dpm_partition_t& partition, sequence_data_t<cluster_tag_t>& data,
          size_t tfbs_length)
{
        // begin with cluster tag 1, since 0 is reserved for the background
        cluster_tag_t k = 1;

        for (dpm_partition_t::const_iterator it = partition.begin();
             it != partition.end(); it++, k++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const seq_index_t& index = static_cast<const seq_index_t&>(**is);
                        for (size_t i = 0; i < tfbs_length; i++) {
                                data[index[0]][index[1]+i] = k;
                        }
                }
        }
}

static void
clean_data(const dpm_partition_t& partition, sequence_data_t<cluster_tag_t>& data,
          size_t tfbs_length)
{
        for (dpm_partition_t::const_iterator it = partition.begin();
             it != partition.end(); it++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const seq_index_t& index = static_cast<const seq_index_t&>(**is);
                        for (size_t i = 0; i < tfbs_length; i++) {
                                data[index[0]][index[1]+i] = 0;
                        }
                }
        }
}

static size_t
distance(const boost::unordered_set<seq_index_t>& indices,
         const sequence_data_t<cluster_tag_t>& a,
         const sequence_data_t<cluster_tag_t>& b)
{
        size_t d = 0;

        for (boost::unordered_set<seq_index_t>::const_iterator it = indices.begin();
             it != indices.end(); it++) {
                for (boost::unordered_set<seq_index_t>::const_iterator is = it;
                     is != indices.end(); is++) {
                        if ((a[*it] == a[*is] && b[*it] != b[*is]) ||
                            (a[*it] != a[*is] && b[*it] == b[*is])) {
                                d += 1;
                        }
                }
        }
        return d;
}

static size_t
distance(const dpm_partition_t& pi_a,
         const dpm_partition_t& pi_b,
         sequence_data_t<cluster_tag_t>& a,
         sequence_data_t<cluster_tag_t>& b,
         size_t tfbs_length)
{
        boost::unordered_set<seq_index_t> indices;

        // initialize auxiliary cluster information
        init_data(pi_a, a, tfbs_length);
        init_data(pi_b, b, tfbs_length);

        // go through subsets of partition a and record all positions
        // that are assigned to a cluster
        for (dpm_partition_t::const_iterator it = pi_a.begin();
             it != pi_a.end(); it++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const seq_index_t& tmp = static_cast<const seq_index_t&>(**is);
                        for (size_t i = 0; i < tfbs_length; i++) {
                                indices.insert(seq_index_t(tmp[0], tmp[1]+i));
                        }
                }
        }
        // go through subsets of partition b and record all positions
        // that are assigned to a cluster
        for (dpm_partition_t::const_iterator it = pi_b.begin();
             it != pi_b.end(); it++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const seq_index_t& tmp = static_cast<const seq_index_t&>(**is);
                        for (size_t i = 0; i < tfbs_length; i++) {
                                indices.insert(seq_index_t(tmp[0], tmp[1]+i));
                        }
                }
        }
        // clean data
        clean_data(pi_a, a, tfbs_length);
        clean_data(pi_b, b, tfbs_length);

        return distance(indices, a, b);
}

static
double mean_loss(double d)
{
        return d*d;
}

static
double median_loss(double d)
{
        return d;
}

static dpm_partition_t
dpm_tfbs_estimate(const dpm_partition_list_t& partitions,
                  const vector<size_t>& sizes, size_t tfbs_length,
                  bool verbose, double (*loss)(double))
{
        // number of partitions
        const size_t n = partitions.size();

        // return if no partitions are given
        if (n == 0) {
                return dpm_partition_t();
        }

        // auxiliary storage
        sequence_data_t<cluster_tag_t> a(sizes, 0);
        sequence_data_t<cluster_tag_t> b(sizes, 0);

        // matrix of distances between every pair of partitions
        matrix<size_t> distances(n, n);
        vector<size_t> sums(n, 0);

        // save value and position of minimum distances
        size_t min = SIZE_MAX, argmin = SIZE_MAX;

        for (size_t i = 0, k = 0; i < n; i++) {
                for (size_t j = i+1; j < n; j++, k++) {
                        if (verbose) {
                                cerr << progress_t(2.0*(k+1)/(double)(n*n-n));
                        }
                        distances[i][j] = distance(partitions[i], partitions[j], a, b, tfbs_length);
                        distances[j][i] = distances[i][j];
                }
        }
        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                        sums[i] += (*loss)(distances[i][j]);
                }
                if (sums[i] < min) {
                        min = sums[i];
                        argmin = i;
                }
        }
        return partitions[argmin];
}

static dpm_partition_t
dpm_tfbs_estimate(const sampling_history_t& history,
                  const vector<size_t>& sizes, size_t tfbs_length,
                  ssize_t take, bool verbose, double (*loss)(double))
{
        /* number of parallel samplers */
        size_t n = history.temperature.size();

        /* if there are no samples, return immediately */
        if (n == 0) return dpm_partition_t();

        /* number of samples for each sampler */
        size_t m = history.temperature[0].size();

        /* list of feasible partitions */
        dpm_partition_list_t partitions;

        /* number of partitions that were discarded */
        size_t discarded = 0;

        for (size_t i = 0; i < m; i++) {
                for (size_t j = 0; j < n; j++) {
                        /* if the temperature is greater than one
                         * discard the sample, since it belongs
                         * to the burnin period */
                        if (history.temperature[j][i] > 1.0) {
                                discarded++;
                                continue;
                        }
                        assert(i*n+j < history.partitions.size());
                        /* save partition */
                        partitions.push_back(history.partitions[i*n+j]);
                }
        }
        /* select a subset of partitions */
        if (take != -1 && (ssize_t)partitions.size() > take) {
                if (verbose) {
                        cerr << boost::format("(Taking only the last %d partitions)") % take
                             << endl;
                }
                partitions = dpm_partition_list_t(partitions.end()-take, partitions.end());
        }
        else {
                if (verbose) {
                        cerr << boost::format("(Discarded first %d partitions)") % discarded
                             << endl;
                }
        }
        /* compute estimate */
        return dpm_tfbs_estimate(partitions, sizes, tfbs_length, verbose, loss);
}

// functions for computing map partitions
// -----------------------------------------------------------------------------

static bool
map_local_optimization(cluster_t& cluster, dpm_tfbs_t& dpm, bool verbose)
{
        double posterior_ref = dpm.posterior();
        double posterior_left;
        double posterior_right;
        stringstream ss;
        size_t size = cluster.size();

        dpm.state().save(cluster.cluster_tag());
        dpm.state().move_left(cluster);
        posterior_left = dpm.posterior();
        dpm.state().restore();

        dpm.state().save(cluster.cluster_tag());
        dpm.state().move_right(cluster);
        posterior_right = dpm.posterior();
        dpm.state().restore();

        if (posterior_left > posterior_ref && posterior_left > posterior_right) {
                dpm.state().move_left(cluster);
                ss << "moved to the left";
        }
        else if (posterior_right > posterior_ref && posterior_right > posterior_left) {
                dpm.state().move_right(cluster);
                ss << "moved to the right";
        }
        else {
                return false;
        }
        if (verbose) {
                cout << boost::format("Cluster %d %s (%d -> %d)")
                        % cluster.cluster_tag()
                        % ss.str()
                        % size
                        % cluster.size()
                     << endl;
        }

        return true;
}

static bool
map_local_optimization(const index_i& index, dpm_tfbs_t& dpm, bool verbose) {
        ////////////////////////////////////////////////////////////////////////
        // check if we can sample this element
        if (!dpm.valid_for_sampling(index)) {
                return false;
        }
        ////////////////////////////////////////////////////////////////////////
        // release the element from its cluster
        cluster_tag_t old_cluster_tag = dpm.state()[index];
        dpm.state().remove(index, old_cluster_tag);
        size_t components = dpm.mixture_components() + dpm.baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        dpm.mixture_weights(index, log_weights, cluster_tags);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        cluster_tag_t new_cluster_tag;
        new_cluster_tag = cluster_tags[select_max_component(components, log_weights)];
        if (verbose && new_cluster_tag != old_cluster_tag) {
                cout << "Moving " << (const seq_index_t&)index
                     << " from cluster " << old_cluster_tag
                     << " to cluster "   << new_cluster_tag
                     << endl;
        }

        ////////////////////////////////////////////////////////////////////////
        dpm.state().add(index, new_cluster_tag);

        return old_cluster_tag != new_cluster_tag;
}

static dpm_partition_t
map_local_optimization(dpm_tfbs_t& dpm, bool verbose) {
        bool optimized;
        /* make some noise */
        if (verbose) {
                cout << "Performing local optimizations..." << endl
                     << "Posterior: " << dpm.posterior()
                     << endl;
        }
        /* loop until no local optimizations are possible */
        do {
                optimized = false;
                /* optimize single positions */
                for (indexer_t::const_iterator it = dpm.data().sampling_begin();
                     it != dpm.data().sampling_end(); it++) {
                        optimized |= map_local_optimization(**it, dpm, verbose);
                }
                /* optimize by moving whole clusters */
                for (mixture_state_t::iterator it = dpm.state().begin();
                     it != dpm.state().end(); it++) {
                        optimized |= map_local_optimization(**it, dpm, verbose);
                }
                if (verbose) {
                        cout << "Posterior: " << dpm.posterior()
                             << endl;
                }
        } while (optimized);
        /* return final partition */
        return dpm.state().partition();
}

static dpm_partition_t
compute_map(const sampling_history_t& history, dpm_tfbs_t& dpm, bool verbose)
{
        /* number of parallel samplers */
        size_t n = history.temperature.size();

        /* if there are no samples, return immediately */
        if (n == 0) return dpm_partition_t();

        /* number of samples for each sampler */
        size_t m = history.temperature[0].size();

        /* number of partitions that were discarded */
        size_t discarded = 0;

        /* max and argmax of posterior samples */
        double max = -numeric_limits<double>::infinity();
        size_t argmax = 0;

        for (size_t i = 0; i < m; i++) {
                for (size_t j = 0; j < n; j++) {
                        /* if the temperature is greater than one
                         * discard the sample, since it belongs
                         * to the burnin period */
                        if (history.temperature[j][i] > 1.0) {
                                discarded++;
                                continue;
                        }
                        assert(i*n+j < history.partitions.size());
                        if (max < history.posterior[j][i]) {
                                max = history.posterior[j][i];
                                argmax = i*n+j;
                        }
                }
        }
        if (verbose) {
                cout << boost::format("(Discarded first %d partitions)") % discarded
                     << endl;
        }
        if (max > -numeric_limits<double>::infinity()) {
                /* set dpm state to best partition */
                dpm.state().set_partition(history.partitions[argmax]);
                /* optimize this partition locally */
                return map_local_optimization(dpm, verbose);
        }
        else {
                return dpm_partition_t();
        }
}

// entry points
// -----------------------------------------------------------------------------

dpm_partition_t
dpm_tfbs_t::map(const sampling_history_t& history, bool verbose) const
{
        /* create a copy of this object */
        dpm_tfbs_t dpm_tfbs(*this);

        if (verbose) {
                cout << "Computing map partition: ";
        }

        return compute_map(history, dpm_tfbs, verbose);
}

dpm_partition_t
dpm_tfbs_t::mean(const sampling_history_t& history, ssize_t take, bool verbose) const
{
        if (verbose) {
                cout << "Computing mean partition: ";
        }
        return dpm_tfbs_estimate(history, data().sizes(), _tfbs_length, take, verbose, &mean_loss);
}

dpm_partition_t
dpm_tfbs_t::median(const sampling_history_t& history, ssize_t take, bool verbose) const
{
        if (verbose) {
                cout << "Computing median partition: ";
        }
        return dpm_tfbs_estimate(history, data().sizes(), _tfbs_length, take, verbose, &median_loss);
}
