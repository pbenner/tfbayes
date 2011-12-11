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

#include <iostream>

#include <dpm-tfbs.hh>
#include <statistics.hh>

using namespace std;

void
DpmTfbs::test_metropolis_hastings() {
        seq_index_t index1(0,53);
        seq_index_t index2(1,26);
        double l1, l2;

        Cluster& cluster1 = _clustermanager.get_free_cluster(_model_tags[0]);
        cluster_tag_t cluster_tag1 = cluster1.cluster_tag();

        cout.precision(10);

        remove(index1, bg_cluster_tag);
        add(index1, cluster_tag1);

        remove(index2, bg_cluster_tag);
        add(index2, cluster_tag1);

        l1 = likelihood();
        cout << _cluster_assignments << endl;
        cout << "likelihood: " << l1 << endl;

        //proposal(cluster1);
        metropolis_hastings();
        l2 = likelihood();
        cout << _cluster_assignments << endl;
        cout << "likelihood: " << l2 << endl;

        cout << "ratio: " << exp(l2-l1) << endl;

        exit(EXIT_SUCCESS);
}

void
DpmTfbs::test_moves() {
        seq_index_t index1(0,0);
        seq_index_t index2(0,10);
        seq_index_t index3(1,10);
        seq_index_t index4(1,20);
        seq_index_t index5(0,13);
        seq_index_t index6(0,22);

        Cluster& cluster1 = _clustermanager.get_free_cluster(_model_tags[0]);
        cluster_tag_t cluster_tag1 = cluster1.cluster_tag();

        remove(index1, bg_cluster_tag);
        add(index1, cluster_tag1);

        remove(index2, bg_cluster_tag);
        add(index2, cluster_tag1);

        remove(index3, bg_cluster_tag);
        add(index3, cluster_tag1);

        Cluster& cluster2 = _clustermanager.get_free_cluster(_model_tags[0]);
        cluster_tag_t cluster_tag2 = cluster2.cluster_tag();

        remove(index4, bg_cluster_tag);
        add(index4, cluster_tag2);

        cout << _cluster_assignments << endl;
        move_left(cluster1);
        cout << _cluster_assignments << endl;

        exit(EXIT_SUCCESS);
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
