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

#define __STDC_LIMIT_MACROS

#include <stdint.h>
#include <cmath> /* abs */

#include <boost/unordered/unordered_set.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/dpm/dpm-sampling-history.hh>
#include <tfbayes/utility/progress.hh>
#include <tfbayes/utility/statistics.hh>

using namespace std;

// utilities for computing means and medians
// -----------------------------------------------------------------------------

static void
init_data(const dpm_partition_t& partition, sequence_data_t<cluster_tag_t>& data)
{
        /* Cluster labeling convention:
         * The kth cluster gets all labels from k*foreground_length + 1 to
         * k*foreground_length + foreground_length. Hence, each site in a
         * cluster gets a unique label. The label 0 is reserved for
         * the background */
        cluster_tag_t k = 0;

        for (dpm_partition_t::const_iterator it = partition.begin();
             it != partition.end(); it++, k++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const index_t& index = is->index();
                        for (size_t i = 0; i < is->length(); i++) {
                                data[index[0]][index[1]+i] = k*is->length()+i+1;
                        }
                }
        }
}

static void
clean_data(const dpm_partition_t& partition, sequence_data_t<cluster_tag_t>& data)
{
        for (dpm_partition_t::const_iterator it = partition.begin();
             it != partition.end(); it++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const index_t& index = is->index();
                        for (size_t i = 0; i < is->length(); i++) {
                                data[index[0]][index[1]+i] = 0;
                        }
                }
        }
}

// naive distance function (for debugging purposes)
// -----------------------------------------------------------------------------

GCC_ATTRIBUTE_UNUSED
static size_t
naive_distance(const dpm_partition_t& pi_a,
               const dpm_partition_t& pi_b,
               sequence_data_t<cluster_tag_t>& a,
               sequence_data_t<cluster_tag_t>& b,
               const dpm_tfbs_t& dpm)
{
        boost::unordered_set<index_t> indices;
        // resulting distance
        size_t d = 0;

        // initialize auxiliary cluster information
        init_data(pi_a, a);
        init_data(pi_b, b);

        for (indexer_t::const_iterator it = dpm.data().begin();
             it != dpm.data().end(); it++) {
                for (indexer_t::const_iterator is = it+1;
                     is != dpm.data().end(); is++) {
                        const index_t& i = *it;
                        const index_t& j = *is;
                        if ((a[i] == a[j] && b[i] != b[j]) ||
                            (a[i] != a[j] && b[i] == b[j])) {
                                d += 1;
                        }
                }
        }

        // clean data
        clean_data(pi_a, a);
        clean_data(pi_b, b);

        return d;
}

// a slightly less naive implementation of a distance function
// bug still quadratic in the number of foreground indices!
// -----------------------------------------------------------------------------

GCC_ATTRIBUTE_UNUSED
static size_t
distance1(const boost::unordered_set<index_t>& indices,
          const sequence_data_t<cluster_tag_t>& a,
          const sequence_data_t<cluster_tag_t>& b,
          size_t bg_size)
{
        size_t d = 0, bg = 0;

        for (boost::unordered_set<index_t>::const_iterator it = indices.begin();
             it != indices.end(); it++) {
                for (boost::unordered_set<index_t>::const_iterator is = it;
                     is != indices.end(); is++) {
                        if ((a[*it] == a[*is] && b[*it] != b[*is]) ||
                            (a[*it] != a[*is] && b[*it] == b[*is])) {
                                d++;
                        }
                }
                if (a[*it] == 0 || b[*it] == 0) {
                        bg++;
                }
        }
        return d + bg*bg_size;
}

// highly optimized distance function (formulas can be found in
// Lawrence Hubert. Comparing Partitions. Journal of
// Classification, 1985)
// -----------------------------------------------------------------------------

static double
distance2(const boost::unordered_set<index_t>& indices,
          const dpm_partition_t& pi_a,
          const dpm_partition_t& pi_b,
          const sequence_data_t<cluster_tag_t>& a,
          const sequence_data_t<cluster_tag_t>& b,
          size_t bg_size)
{
        // boost sparse matrices
        using namespace boost::numeric::ublas;

        // resulting distance
        double result = 0;

        // number of clusters
        size_t la = 1; // start at one for the
        size_t lb = 1; // background model
        for (dpm_partition_t::const_iterator it = pi_a.begin();
             it != pi_a.end(); it++) {
                la += it->begin()->length();
        }
        for (dpm_partition_t::const_iterator it = pi_b.begin();
             it != pi_b.end(); it++) {
                lb += it->begin()->length();
        }

        // contingency table
        mapped_matrix<double> m(la, lb);
        std::vector<double> ma(la, 0.0);
        std::vector<double> mb(lb, 0.0);

        for (boost::unordered_set<index_t>::const_iterator it = indices.begin();
             it != indices.end(); it++) {
                m(a[*it], b[*it]) += 1.0;
                ma[a[*it]] += 1.0;
                mb[b[*it]] += 1.0;
        }

        // compute distance from contingency table
        for (size_t i = 0; i < la; i++) {
                result += 1.0/2.0*ma[i]*ma[i];
        }
        for (size_t j = 0; j < lb; j++) {
                result += 1.0/2.0*mb[j]*mb[j];
        }
        for (mapped_matrix<double>::const_iterator1 it = m.begin1();
             it != m.end1(); it++) {
                for (mapped_matrix<double>::const_iterator2 is = it.begin();
                     is != it.end(); is++) {
                        result -= (*is)*(*is);
                }
        }
        return result + (ma[0]+mb[0])*bg_size;
}

static double
distance(const dpm_partition_t& pi_a,
         const dpm_partition_t& pi_b,
         sequence_data_t<cluster_tag_t>& a,
         sequence_data_t<cluster_tag_t>& b,
         const dpm_tfbs_t& dpm)
{
        // set of indices that are actually used in clusters
        boost::unordered_set<index_t> indices;

        // initialize auxiliary cluster information
        init_data(pi_a, a);
        init_data(pi_b, b);

        // go through subsets of partition a and record all positions
        // that are assigned to a cluster
        for (dpm_partition_t::const_iterator it = pi_a.begin();
             it != pi_a.end(); it++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const index_t& tmp = is->index();
                        for (size_t i = 0; i < is->length(); i++) {
                                indices.insert(index_t(tmp[0], tmp[1]+i));
                        }
                }
        }
        // go through subsets of partition b and record all positions
        // that are assigned to a cluster
        for (dpm_partition_t::const_iterator it = pi_b.begin();
             it != pi_b.end(); it++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const index_t& tmp = is->index();
                        for (size_t i = 0; i < is->length(); i++) {
                                indices.insert(index_t(tmp[0], tmp[1]+i));
                        }
                }
        }
        // size of the background cluster
        size_t bg_size = dpm.data().elements() - indices.size();

        // compute the distance
//        size_t d = distance1(indices, a, b, bg_size);
        double d = distance2(indices, pi_a, pi_b, a, b, bg_size);

        // clean data
        clean_data(pi_a, a);
        clean_data(pi_b, b);

        return d;
}

static double
distance(const dpm_partition_list_t& partitions,
         const dpm_partition_t& pi,
         sequence_data_t<cluster_tag_t>& a,
         sequence_data_t<cluster_tag_t>& b,
         const dpm_tfbs_t& dpm)
{
        double d = 0.0;

        for (dpm_partition_list_t::const_iterator it = partitions.begin();
             it != partitions.end(); it++) {
                d += distance(*it, pi, a, b, dpm);
        }
        return d;
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
                  double (*loss)(double),
                  const dpm_tfbs_t& dpm,
                  bool verbose)
{
        // number of partitions
        const size_t n = partitions.size();

        // return if no partitions are given
        if (n == 0) {
                return dpm_partition_t();
        }

        // auxiliary storage
        sequence_data_t<cluster_tag_t> a(dpm.data().sizes(), 0);
        sequence_data_t<cluster_tag_t> b(dpm.data().sizes(), 0);

        // matrix of distances between every pair of partitions
        matrix<size_t> distances(n, n);
        vector<size_t> sums(n, 0);

        // save value and position of minimum distances
        size_t min = SIZE_MAX, argmin = SIZE_MAX;

        for (size_t i = 0, k = 0; i < n; i++) {
                for (size_t j = i+1; j < n; j++, k++) {
                        if (verbose && ((k+1) % 100 == 0 || 2*(k+1) == n*n-n)) {
                                cerr << progress_t(2.0*(k+1)/(double)(n*n-n));
                        }
                        distances[i][j] = distance(partitions[i], partitions[j], a, b, dpm);
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

static bool
dpm_tfbs_optimize_estimate(const dpm_partition_list_t& partitions,
                           dpm_partition_t& estimate,
                           double (*loss)(double),
                           sequence_data_t<cluster_tag_t>& a,
                           sequence_data_t<cluster_tag_t>& b,
                           const dpm_tfbs_t& dpm,
                           bool verbose)
{
        // keep track of losses
        double old_loss = (*loss)(distance(partitions, estimate, a, b, dpm));
        double new_loss = old_loss;

        // did we optimize so far?
        bool optimized = false;

        // run optimization
        for (dpm_partition_t::iterator it = estimate.begin();
             it != estimate.end(); it++) {
                if (verbose) {
                        cerr << progress_t((it - estimate.begin() + 1)/(double)(estimate.end()-estimate.begin()));
                }

                // reference to current subset
                dpm_subset_t& subset = *it;

                // current size of the subset
                size_t subset_size = subset.size();

                // set of indices that need to be optimized
                boost::unordered_set<range_t> range_set;

                // populate index list
                for (dpm_subset_t::iterator is = subset.begin();
                     is != subset.end(); is++) {
                        range_set.insert(*is);
                }

                for (boost::unordered_set<range_t>::const_iterator is = range_set.begin();
                     is != range_set.end(); is++) {
                        // save index pointer
                        const range_t range(*is);
                        // erase index from subset
                        subset.erase(range);
                        // compute new loss
                        new_loss = (*loss)(distance(partitions, estimate, a, b, dpm));
                        if (new_loss < old_loss && abs(new_loss - old_loss) > 1e-5) {
                                // new estimate seems better
                                old_loss = new_loss;
                                optimized = true;
                        }
                        else {
                                // new estimate is worse
                                new_loss = old_loss;
                                // reverse changes
                                subset.insert(range);
                        }
                }
                // check if subset is empty
                if (subset.size() == 0) {
                        it = estimate.erase(it);
                        if (verbose) {
                                cout << boost::format("\rRemoved cluster with %d element(s).")
                                        % subset_size
                                     << "\033[K" // clear till end of line
                                     << endl;
                        }
                }
        }

        return optimized;
}

static dpm_partition_t
dpm_tfbs_optimize_estimate(const dpm_partition_list_t& partitions,
                           const dpm_partition_t& estimate,
                           double (*loss)(double),
                           const dpm_tfbs_t& dpm,
                           bool verbose)
{
        dpm_partition_t current_estimate(estimate);

        // return if no partitions are given
        if (partitions.size() == 0) {
                return dpm_partition_t();
        }

        // auxiliary storage
        sequence_data_t<cluster_tag_t> a(dpm.data().sizes(), 0);
        sequence_data_t<cluster_tag_t> b(dpm.data().sizes(), 0);

        // record if distance was minimized
        bool optimized;

        /* make some noise */
        if (verbose) {
                cout << "Performing local optimizations..." << endl;
        }
        /* loop until no local optimizations are possible */
        do {
                /* optimize by removing elements from clusters */
                optimized = dpm_tfbs_optimize_estimate(partitions, current_estimate, loss, a, b, dpm, verbose);

        } while (optimized);

        return current_estimate;
}

static dpm_partition_t
dpm_tfbs_estimate(const sampling_history_t& history,
                  ssize_t take,
                  double (*loss)(double),
                  dpm_tfbs_t& dpm,
                  bool verbose)
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
        dpm_partition_t estimate = dpm_tfbs_estimate(partitions, loss, dpm, verbose);
        /* optimize estimate */
        estimate = dpm_tfbs_optimize_estimate(partitions, estimate, loss, dpm, verbose);

        return estimate;
}

// functions for computing map partitions
// -----------------------------------------------------------------------------

static bool
map_local_optimization_move(cluster_t& cluster, cluster_tag_t bg_cluster_tag, dpm_tfbs_t& dpm, bool verbose)
{
        double posterior_ref = dpm.posterior();
        double posterior_left;
        double posterior_right;
        stringstream ss;
        size_t size = cluster.size();

        dpm.state().save(cluster.cluster_tag(), bg_cluster_tag);
        dpm.state().move_left(cluster, bg_cluster_tag);
        posterior_left = dpm.posterior();
        dpm.state().restore();

        dpm.state().save(cluster.cluster_tag(), bg_cluster_tag);
        dpm.state().move_right(cluster, bg_cluster_tag);
        posterior_right = dpm.posterior();
        dpm.state().restore();

        if (posterior_left > posterior_ref && posterior_left > posterior_right) {
                dpm.state().move_left(cluster, bg_cluster_tag);
                ss << "moved to the left";
        }
        else if (posterior_right > posterior_ref && posterior_right > posterior_left) {
                dpm.state().move_right(cluster, bg_cluster_tag);
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
map_local_optimization_block(cluster_t& cluster, dpm_tfbs_t& dpm, bool verbose)
{
        vector<range_t> range_set;
        cluster_tag_t old_cluster_tag = cluster.cluster_tag();
        bool result = false;

        ////////////////////////////////////////////////////////////////////////
        // fill range_set
        for (cluster_t::iterator it = cluster.begin(); it != cluster.end(); it++)
        {
                const range_t& range = *it;

                range_set.push_back(range);
        }
        ////////////////////////////////////////////////////////////////////////
        // release all elemente from the cluster
        for (vector<range_t>::iterator it = range_set.begin(); it != range_set.end(); it++)
        {
                const range_t& range = *it;

                dpm.state().remove(range);
        }
        ////////////////////////////////////////////////////////////////////////
        // obtian the mixture probabilities
        size_t components = dpm.mixture_components() + dpm.baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        dpm.mixture_weights(range_set, log_weights, cluster_tags, 1.0, true);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        cluster_tag_t new_cluster_tag = cluster_tags[select_max_component(components, log_weights)];

        ////////////////////////////////////////////////////////////////////////
        // print some information to stdout
        flockfile(stdout);
        if (dpm.state()[new_cluster_tag].size() != 0) {
                cout << boost::format("Cluster %d merged with cluster %d (%d + %d)")
                        % old_cluster_tag % new_cluster_tag
                        % dpm.state()[new_cluster_tag].size() % range_set.size()
                     << endl;
                result = true;
        }
        fflush(stdout);
        funlockfile(stdout);

        ////////////////////////////////////////////////////////////////////////
        // move all elements to the new cluster
        for (size_t i = 0; i < range_set.size(); i++) {
                dpm.state().add(range_set[i], new_cluster_tag);
        }
        return result;
}

static bool
map_local_optimization(cluster_t& cluster, dpm_tfbs_t& dpm, bool verbose)
{
        bool result = false;

        for (dpm_tfbs_state_t::bg_cluster_tags_t::const_iterator it = dpm.state().bg_cluster_tags.begin();
             it != dpm.state().bg_cluster_tags.end(); it++) {
                result |= map_local_optimization_move(cluster, *it, dpm, verbose);
        }
        result |= map_local_optimization_block(cluster, dpm, verbose);

        return result;
}

static bool
map_local_optimization(const index_t& index, dpm_tfbs_t& dpm, bool verbose) {
        size_t length;
        if (dpm.state().get_free_range(index, length)) {
                return false;
        }
        range_t range1(index, length, false);
        range_t range2(index, length, true );
        ////////////////////////////////////////////////////////////////////////
        // first release the element from its cluster
        cluster_tag_t old_cluster_tag = dpm.state()[index];
        cluster_tag_t new_cluster_tag;
        dpm.state().remove(range1);
        ////////////////////////////////////////////////////////////////////////
        size_t components = dpm.mixture_components() + dpm.baseline_components();
        double log_weights1[components];
        double log_weights2[components];
        cluster_tag_t cluster_tags1[components];
        cluster_tag_t cluster_tags2[components];
        ////////////////////////////////////////////////////////////////////////
        // compute weights
        dpm.mixture_weights(range1, log_weights1, cluster_tags1, 1.0, true);
        dpm.mixture_weights(range2, log_weights2, cluster_tags2, 1.0, false, log_weights1[components-1]);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        std::pair<size_t, size_t> result = select_max_component2(components, log_weights1, log_weights2);
        ////////////////////////////////////////////////////////////////////////
        if (result.first == 1) {
                new_cluster_tag = cluster_tags1[result.second];
                dpm.state().add(range1, new_cluster_tag);
        }
        else {
                new_cluster_tag = cluster_tags2[result.second];
                dpm.state().add(range2, new_cluster_tag);
        }
        if (verbose && new_cluster_tag != old_cluster_tag) {
                cout << "Moving " << (const index_t&)index
                     << " from cluster " << old_cluster_tag
                     << " to cluster "   << new_cluster_tag
                     << endl;
        }
        return old_cluster_tag != new_cluster_tag;
}

static dpm_partition_t
map_local_optimization(dpm_tfbs_t& dpm, bool verbose) {
        bool optimized;
        double old_posterior;
        double new_posterior = dpm.posterior();
        /* make some noise */
        if (verbose) {
                cout << "Performing local optimizations..." << endl
                     << "Posterior: " << new_posterior
                     << endl;
        }
        /* loop until no local optimizations are possible */
        do {
                optimized = false;
                /* optimize single positions */
                for (indexer_t::const_iterator it = dpm.data().sampling_begin();
                     it != dpm.data().sampling_end(); it++) {
                        optimized |= map_local_optimization(*it, dpm, verbose);
                }
                // since clusters are modified it is not possible to simply
                // loop through the list of clusters, we need to be a bit more
                // careful here!
                vector<cluster_tag_t> used_clusters;
                for (mixture_state_t::const_iterator it = dpm.state().begin(); it != dpm.state().end(); it++) {
                        cluster_t& cluster = **it;
                        if (!dpm.state().is_background(cluster)) {
                                used_clusters.push_back(cluster.cluster_tag());
                        }
                }
                // go through the list of used clusters and if they are still
                // used then generate a block sample
                for (vector<cluster_tag_t>::const_iterator it = used_clusters.begin(); it != used_clusters.end(); it++) {
                        cluster_t& cluster = dpm.state()[*it];
                        if (cluster.size() != 0) {
                                optimized |= map_local_optimization(cluster, dpm, verbose);
                        }
                }
                old_posterior = new_posterior;
                new_posterior = dpm.posterior();

                if (verbose) {
                        cout << "Posterior: " << new_posterior
                             << " (increment: " << abs(old_posterior - new_posterior) << ")"
                             << endl;
                }
        } while (optimized && abs(old_posterior - new_posterior) > 1e-4);
        /* return final partition */
        return dpm.state().partition();
}

static dpm_partition_t
compute_map(const sampling_history_t& history, dpm_tfbs_t& dpm, bool optimize, bool verbose)
{
        /* number of parallel samplers */
        size_t n = history.temperature.size();

        /* if there are no samples, return immediately */
        if (n == 0) return dpm_partition_t();

        /* number of samples for each sampler */
        size_t m = history.temperature[0].size();

        assert(m*n == history.partitions.size());

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
                cout << boost::format("Posterior: %f") % dpm.posterior()
                     << endl;
                if (optimize) {
                        /* optimize this partition locally */
                        return map_local_optimization(dpm, verbose);
                }
                else {
                        return history.partitions[argmax];
                }
        }
        else {
                return dpm_partition_t();
        }
}

// Entry points
// -----------------------------------------------------------------------------

dpm_partition_t
dpm_tfbs_t::map(const sampling_history_t& history, bool optimize, bool verbose) const
{
        /* create a copy of this object */
        dpm_tfbs_t dpm(*this);

        if (verbose) {
                cout << "Computing map partition: ";
        }
        return compute_map(history, dpm, optimize, verbose);
}

dpm_partition_t
dpm_tfbs_t::mean(const sampling_history_t& history, ssize_t take, bool verbose) const
{
        /* create a copy of this object */
        dpm_tfbs_t dpm(*this);

        if (verbose) {
                cout << "Computing mean partition: ";
        }
        return dpm_tfbs_estimate(history, take, &mean_loss, dpm, verbose);
}

dpm_partition_t
dpm_tfbs_t::median(const sampling_history_t& history, ssize_t take, bool verbose) const
{
        /* create a copy of this object */
        dpm_tfbs_t dpm(*this);

        if (verbose) {
                cout << "Computing median partition: ";
        }
        return dpm_tfbs_estimate(history, take, &median_loss, dpm, verbose);
}
