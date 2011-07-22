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

#ifndef CLUSTERMANAGER_HH
#define CLUSTERMANAGER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <vector>

#include <cluster.hh>
#include <datatypes.hh>
#include <distribution.hh>

// ClusterManager
////////////////////////////////////////////////////////////////////////////////

// The cluster manager has a list of clusters and keeps track which of
// them are used. If needed it allocates more free clusters.

class ClusterManager : public Observer<cluster_event_t> {
public:
         ClusterManager(Distribution* distribution, data_t<cluster_tag_t>& cluster_assignments);
         ClusterManager(const ClusterManager& cm);
        ~ClusterManager();

        // iterators
        ////////////////////////////////////////////////////////////////////////
        typedef std::list<Cluster*>::iterator iterator;
        typedef std::list<Cluster*>::const_iterator const_iterator;

        iterator begin() { return used_clusters.begin(); }
        iterator end()   { return used_clusters.end();   }

        const_iterator begin() const { return used_clusters.begin(); }
        const_iterator end()   const { return used_clusters.end();   }

        typedef std::vector<Cluster*>::iterator iterator_all;
        typedef std::vector<Cluster*>::const_iterator const_iterator_all;

        iterator_all begin_all() { return clusters.begin(); }
        iterator_all end_all()   { return clusters.end();   }

        const_iterator_all begin_all() const { return clusters.begin(); }
        const_iterator_all end_all()   const { return clusters.end();   }

        // operators
        ////////////////////////////////////////////////////////////////////////
        __inline__       Cluster& operator[](cluster_tag_t c)            { return *clusters[c];               }
        __inline__ const Cluster& operator[](cluster_tag_t c)      const { return *clusters[c];               }
        __inline__  cluster_tag_t operator[](const index_t& index) const { return cluster_assignments[index]; }

        friend std::ostream& operator<< (std::ostream& o, const ClusterManager& cm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        void update(Observed<cluster_event_t>* cluster, cluster_event_t event);
        void update(Observed<cluster_event_t>* cluster, cluster_event_t event, const range_t& range);
        cluster_tag_t add_cluster();
        cluster_tag_t add_cluster(Distribution* distribution);
        Cluster& get_free_cluster();
        __inline__ size_t size() const { return used_clusters_size; };

        // keep track of cluster assignments
        void record_cluster_assignment(const range_t& range, cluster_tag_t tag);

private:
        std::vector<Cluster*> clusters;
        std::list<Cluster*> used_clusters;
        std::list<Cluster*> free_clusters;

        // list.size() is inefficient, so keep track of
        // the number of elements
        size_t used_clusters_size;
        size_t free_clusters_size;

        // default distribution for the dirichlet process
        Distribution* default_distribution;

        // assignments to clusters
        data_t<cluster_tag_t>& cluster_assignments;
};

#endif /* CLUSTERMANAGER_HH */
