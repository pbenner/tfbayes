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

ClusterManager::ClusterManager()
{
}

ClusterManager::ClusterManager(Distribution* distribution, data_t<cluster_tag_t>& cluster_assignments)
        : used_clusters_size(0), free_clusters_size(0),
          default_distribution(distribution),
          cluster_assignments(cluster_assignments)
{
}

ClusterManager::~ClusterManager() {
        for (ClusterManager::iterator_all it = begin_all();
             it != end_all(); it++) {
                delete(*it);
        }
        delete(default_distribution);
}

ClusterManager::ClusterManager(const ClusterManager& cm)
        : used_clusters_size(cm.used_clusters_size),
          free_clusters_size(cm.free_clusters_size),
          default_distribution(cm.default_distribution->clone())
{
        for (size_t i = 0; i < cm.clusters.size(); i++) {
                Cluster* c = new Cluster(*cm.clusters[i]);
                c->set_observer(this);
                clusters.push_back(c);
                if (c->size() == 0 && c->destructible()) {
                        free_clusters.push_back(c);
                }
                else {
                        used_clusters.push_back(c);
                }
        }
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
        Cluster* c = new Cluster(distribution, tag, this, false);
        clusters.push_back(c);
        used_clusters.push_back(clusters.back());
        used_clusters_size++;

        return tag;
}

// add cluster with default distribution
cluster_tag_t
ClusterManager::add_cluster()
{
        cluster_tag_t tag = clusters.size();

        Cluster* c = new Cluster(default_distribution->clone(), tag, this);
        clusters.push_back(c);
        // this cluster is empty, so place it into the list
        // of free clusters
        free_clusters.push_back(clusters.back());
        free_clusters_size++;

        return tag;
}

Cluster&
ClusterManager::get_free_cluster() {
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
                if (cluster->destructible()) {
                        used_clusters.remove(cluster);
                        free_clusters.push_back(cluster);
                        used_clusters_size--;
                        free_clusters_size++;
                        break;
                }
        case cluster_event_nonempty:
                if (cluster->destructible()) {
                        used_clusters.push_back(cluster);
                        free_clusters.remove(cluster);
                        used_clusters_size++;
                        free_clusters_size--;
                        break;
                }
        default:
                break;
        }
}

void
ClusterManager::update(Observed<cluster_event_t>* observed, cluster_event_t event, const range_t& range)
{
        Cluster* cluster = (Cluster*)observed;
        iterator_t<cluster_tag_t> iterator = cluster_assignments[range];

        switch (event) {
        case cluster_event_add_word:
                do {
                        *iterator = cluster->tag();
                } while(iterator++);
                break;
        case cluster_event_remove_word:
                do {
                        *iterator = -1;
                } while(iterator++);
                break;
        default:
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

cluster_tag_t
ClusterManager::operator[](const index_t& index) const {
        return cluster_assignments[index];
}

size_t
ClusterManager::size() const {
        return used_clusters_size;
}

// ostream&
// operator<< (ostream& o, const ClusterManager& cm) {
//         for (size_t i = 0; i < cm.cluster_assignments.size(); i++) {
//                 for (size_t j = 0; j < cm.cluster_assignments[i].size(); j++) {
//                         o << cm.cluster_assignments[i][j] << " ";
//                 }
//                 o << endl;
//         }

//         return o;
// }
