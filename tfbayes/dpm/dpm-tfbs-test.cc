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

#include <iostream>
#include <limits>
#include <sstream>
#include <sys/time.h>

#include <boost/random/mersenne_twister.hpp>

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/utility/statistics.hh>

using namespace std;

void
dpm_tfbs_t::test_posterior(cluster_t& cluster, const range_t& range)
{
        double tmp1 = 0;
        double tmp2 = 0;

        if (state().is_background(cluster)) {
                tmp1 = background_mixture_weight(range, cluster);
        }
        else {
                if (cluster.model().id().length > range.length()) {
                        return;
                }
                tmp1 = foreground_mixture_weight(range, cluster);
        }
        state().add(range, cluster.cluster_tag());
        tmp2 += posterior();
        state().remove(range);
        tmp2 -= posterior();
        if (abs(tmp1-tmp2) > 0.001) {
                cerr << "ERROR: posterior probability is inconsistent!" << endl;
                cerr << "-> tmp1: " << tmp1 << endl;
                cerr << "-> tmp2: " << tmp2 << endl;
                exit(EXIT_FAILURE);
        }
}

void
dpm_tfbs_t::test_moves() {
        range_t range1(index_t(0, 0), 10, true);
        range_t range2(index_t(0,10), 10, true);
        range_t range3(index_t(1,10), 10, true);
        range_t range4(index_t(1,20), 10, true);
        range_t range5(index_t(0,13), 10, true);
        range_t range6(index_t(0,22), 10, true);

        cluster_t& cluster1 = m_state.get_free_cluster(m_baseline_tags[0]);
        cluster_tag_t cluster_tag1 = cluster1.cluster_tag();

        m_state.remove(range1);
        m_state.add(range1, cluster_tag1);

        m_state.remove(range2);
        m_state.add(range2, cluster_tag1);

        m_state.remove(range3);
        m_state.add(range3, cluster_tag1);

        cluster_t& cluster2 = m_state.get_free_cluster(m_baseline_tags[0]);
        cluster_tag_t cluster_tag2 = cluster2.cluster_tag();

        m_state.remove(range4);
        m_state.add(range4, cluster_tag2);

        cout << m_state.cluster_assignments() << endl;
        m_state.move_left(cluster1, m_state.bg_cluster_tags[0]);
        cout << m_state.cluster_assignments() << endl;

        exit(EXIT_SUCCESS);
}

void
dpm_tfbs_t::test_background() {
        range_t range1(index_t(0, 0), 10, true);
        range_t range2(index_t(0,10), 10, true);
        range_t range3(index_t(1,10), 10, true);
        range_t range4(index_t(1,20), 10, true);
        range_t range5(index_t(0,13), 10, true);
        range_t range6(index_t(0,22), 10, true);

        cluster_t& cluster1 = m_state.get_free_cluster(m_baseline_tags[0]);
        cluster_tag_t cluster_tag1 = cluster1.cluster_tag();
        cout << "Adding range1 to cluster:" << cluster_tag1 << endl;
        m_state.remove(range1);
        m_state.add(range1, cluster_tag1);
        cout << "Adding range6 to cluster:" << cluster_tag1 << endl;
        m_state.remove(range6);
        m_state.add(range6, cluster_tag1);
        cout << endl;

        cout << m_state.cluster_assignments() << endl;

        cluster_t& cluster2 = m_state.get_free_cluster(m_baseline_tags[0]);
        cluster_tag_t cluster_tag2 = cluster2.cluster_tag();

        cout << "Adding range2 to cluster:" << cluster_tag2 << endl;
        m_state.remove(range2);
        m_state.add(range2, cluster_tag2);
        cout << m_state.cluster_assignments();
        m_state.remove(range2);
        m_state.add(range2, m_state.bg_cluster_tags[0]);
        cout << "Removing range2 from cluster:" << cluster_tag2 << endl << endl;

        cout << "Adding range3 to cluster:" << cluster_tag2 << endl;
        m_state.remove(range3);
        m_state.add(range3, cluster_tag2);
        cout << m_state.cluster_assignments();
        m_state.remove(range3);
        m_state.add(range3, m_state.bg_cluster_tags[0]);
        cout << "Removing range3 from cluster:" << cluster_tag2 << endl << endl;

        cout << "Adding range4 to cluster:" << cluster_tag2 << endl;
        m_state.remove(range4);
        m_state.add(range4, cluster_tag2);
        cout << m_state.cluster_assignments();
        m_state.remove(range4);
        m_state.add(range4, m_state.bg_cluster_tags[0]);
        cout << "Removing range4 from cluster:" << cluster_tag2 << endl << endl;

        cout << "Adding range5 to cluster:" << cluster_tag2 << endl;
        m_state.remove(range5);
        m_state.add(range5, cluster_tag2);
        cout << m_state.cluster_assignments();
        m_state.remove(range5);
        m_state.add(range5, m_state.bg_cluster_tags[0]);
        cout << "Removing range5 from cluster:" << cluster_tag2 << endl;

        exit(EXIT_SUCCESS);
}

#include <tfbayes/utility/logarithmetic.hh>

void normalize_weights(size_t components, double *log_weights)
{
        double sum = -numeric_limits<double>::infinity();

        for (size_t i = 0; i < components; i++) {
                sum = logadd(sum, log_weights[i]);
        }

        for (size_t i = 0; i < components; i++) {
                log_weights[i] -= sum;
        }
}

void
dpm_tfbs_t::test() {
        boost::random::mt19937 gen;
        struct timeval tv;
        gettimeofday(&tv, NULL);
        gen.seed(tv.tv_sec*tv.tv_usec);

        range_t range1(index_t(0,0), 10, true);
        range_t range2(index_t(1,0), 10, true);
        range_t range3(index_t(2,0), 10, true);
        range_t range4(index_t(3,0), 10, true);

        cluster_t& cluster1 = m_state.get_free_cluster(m_baseline_tags[0]);
        cluster_tag_t cluster_tag1 = cluster1.cluster_tag();
        cout << "Adding range1 to cluster:" << cluster_tag1 << endl;
        m_state.add(range1, cluster_tag1);
        cluster_t& cluster2 = m_state.get_free_cluster(m_baseline_tags[0]);
        cluster_tag_t cluster_tag2 = cluster2.cluster_tag();
        cout << "Adding range2 to cluster:" << cluster_tag2 << endl;
        m_state.add(range2, cluster_tag2);

        cout << "Components: " << mixture_components() << " + " << baseline_components() << endl;
        size_t components = mixture_components() + baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        cluster_tag_t new_cluster_tag;

        cout << "Sampling range3" << endl;
        mixture_weights(range3, log_weights, cluster_tags);
        normalize_weights(components, log_weights);
        for (size_t i = 0; i < components; i++) {
                cout << "weight " << i << ": " << exp(log_weights[i]) << endl;
        }
        for (size_t i = 0; i < 100; i++) {
                new_cluster_tag = cluster_tags[select_component(components, log_weights, gen)];
                cout << "selected cluster " << new_cluster_tag << endl;
        }

        cout << "Sampling range4" << endl;
        mixture_weights(range4, log_weights, cluster_tags);
        for (size_t i = 0; i < components; i++) {
                cout << "weight " << i << ": " << log_weights[i] << endl;
        }
        for (size_t i = 0; i < 100; i++) {
                new_cluster_tag = cluster_tags[select_component(components, log_weights, gen)];
                cout << "selected cluster " << new_cluster_tag << endl;
        }

        exit(EXIT_SUCCESS);
}
