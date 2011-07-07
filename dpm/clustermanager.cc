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
#include <cstdlib>
#include <ctime>
#include <list>

#include <clustermanager.hh>

using namespace std;

ClusterManager::ClusterManager(Distribution* distribution)
        : used_clusters_size(0), free_clusters_size(0),
          default_distribution(distribution)
{
}

ClusterManager::~ClusterManager() {
        for (ClusterManager::iterator_all it = begin_all();
             it != end_all(); it++) {
                delete(*it);
        }
        delete(default_distribution);
}

// add a cluster with specific distribution which will
// not be destructible (for other components of the mixture model)
cluster_tag_t
ClusterManager::add_cluster(Distribution* distribution)
{
        cluster_tag_t tag = clusters.size();

        // since this cluster is a fixed component of the
        // process, it should never be put into the list
        // of free clusters, so we do not need to observe it
        // and we can simply push it to the list of used
        // clusters
        Cluster* c = new Cluster(distribution, tag, false);
        clusters.push_back(c);
        used_clusters.push_back(clusters.back());
        used_clusters_size++;

        cout << *(clusters[0]) << endl;

        return tag;
}

// add cluster with default distribution
cluster_tag_t
ClusterManager::add_cluster()
{
        cluster_tag_t tag = clusters.size();

        Cluster* c = new Cluster(default_distribution->clone(), tag, true, this);
        clusters.push_back(c);
        // this cluster is empty, so place it into the list
        // of free clusters
        free_clusters.push_back(clusters.back());
        free_clusters_size++;

        return tag;
}

Cluster&
ClusterManager::get_free_cluster() {
        // TODO: don't call free_clusters.size()
        if (free_clusters_size == 0) {
                // create new cluster
                add_cluster();
        }
        return *free_clusters.front();
}

void
ClusterManager::update(Observed<cluster_event_t>* observed, cluster_event_t event)
{
        Cluster* cluster = (Cluster*)observed;

        switch (event) {
        case cluster_event_empty:
                used_clusters.remove(cluster);
                free_clusters.push_back(cluster);
                used_clusters_size--;
                free_clusters_size++;
                break;
        case cluster_event_nonempty:
                used_clusters.push_back(cluster);
                free_clusters.remove(cluster);
                used_clusters_size++;
                free_clusters_size--;
                break;
        }
}

Cluster&
ClusterManager::operator[](cluster_tag_t c) {
        return *clusters[c];
}

const Cluster&
ClusterManager::operator[](cluster_tag_t c) const {
        return *clusters[c];
}

// ostream&
// operator<< (ostream& o, ClusterManager const& clusters)
// {
//         for (ClusterManager::const_iterator it1 = clusters.begin();
//              it1 != clusters.end(); it1++) {
//                 o << "ClusterManager: " << (*it1)->tag << endl;
//                 // for (ClusterManager::elements_t::const_iterator it2 = (*it1)->elements.begin();
//                 //      it2 != (*it1)->elements.end(); it2++) {
//                 //         if (it2 != (*it1)->elements.begin()) {
//                 //                 o << ", ";
//                 //         }
//                 //         o << **it2;
//                 // }
//                 o << endl << endl;
//         }

//         o << "Used clusters: ";
//         for (list<ClusterManager::cluster*>::const_iterator it = clusters.used_clusters.begin();
//              it != clusters.used_clusters.end(); it++) {
//                 if (it != clusters.used_clusters.begin()) {
//                         o << ", ";
//                 }
//                 o << (*it)->tag;
//         }
//         o << endl;

//         o << "Free clusters: ";
//         for (list<ClusterManager::cluster*>::const_iterator it = clusters.free_clusters.begin();
//              it != clusters.free_clusters.end(); it++) {
//                 if (it != clusters.free_clusters.begin()) {
//                         o << ", ";
//                 }
//                 o << (*it)->tag;
//         }
//         o << endl;

//         return o;
// }
