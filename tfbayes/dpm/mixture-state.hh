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

#ifndef MIXTURE_STATE_HH
#define MIXTURE_STATE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <map>
#include <vector>

#include <tfbayes/dpm/clonable.hh>
#include <tfbayes/dpm/cluster.hh>
#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/dpm/datatypes.hh>
#include <tfbayes/dpm/dpm-partition.hh>

// mixture_state_t
////////////////////////////////////////////////////////////////////////////////

// The cluster manager has a list of clusters and keeps track which of
// them are used. If needed it allocates more free clusters.

class mixture_state_t : public Observer<cluster_event_t>, public clonable {
public:
         mixture_state_t(const data_t<cluster_tag_t>& cluster_assignments);
         mixture_state_t(const mixture_state_t& cm);
        ~mixture_state_t();

        mixture_state_t* clone() const;

        friend void swap(mixture_state_t& first, mixture_state_t& second);

        // iterators
        ////////////////////////////////////////////////////////////////////////
        typedef std::list<cluster_t*>::iterator iterator;
        typedef std::list<cluster_t*>::const_iterator const_iterator;

        iterator begin() { return used_clusters.begin(); }
        iterator end()   { return used_clusters.end();   }

        const_iterator begin() const { return used_clusters.begin(); }
        const_iterator end()   const { return used_clusters.end();   }

        typedef std::vector<cluster_t*>::iterator iterator_all;
        typedef std::vector<cluster_t*>::const_iterator const_iterator_all;

        // operators
        ////////////////////////////////////////////////////////////////////////
        __inline__       cluster_t& operator[](cluster_tag_t c)            { return *clusters[c]; }
        __inline__ const cluster_t& operator[](cluster_tag_t c)      const { return *clusters[c]; }
        __inline__   cluster_tag_t  operator[](const index_i& index) const { return  cluster_assignments()[index]; }

        virtual mixture_state_t& operator=(const mixture_state_t& mixture_state);

        // methods
        ////////////////////////////////////////////////////////////////////////
        void update(Observed<cluster_event_t>* cluster, cluster_event_t event);
        void update(Observed<cluster_event_t>* cluster, cluster_event_t event, const range_t& range);
        void add_baseline_model(component_model_t* distribution, const baseline_tag_t& baseline_tag);
        cluster_tag_t add_cluster(baseline_tag_t baseline_tag);
        cluster_tag_t add_cluster(component_model_t* distribution);
        cluster_t& get_free_cluster(baseline_tag_t baseline_tag);
        __inline__ size_t size() const { return used_clusters_size; }
        // access to cluster assignments, the data type might be
        // different in child classes, so make this virtual
        virtual data_t<cluster_tag_t>& cluster_assignments() { return *_cluster_assignments; }
        virtual const data_t<cluster_tag_t>& cluster_assignments() const { return *_cluster_assignments; }
        // also provide a function to generate a partition from
        // cluster assignments
        virtual dpm_partition_t partition() const;

protected:
        std::vector<cluster_t*> clusters;
        std::list<cluster_t*> used_clusters;
        std::list<cluster_t*> free_clusters;

        // list.size() is inefficient, so keep track of
        // the number of elements
        size_t used_clusters_size;
        size_t free_clusters_size;

        // distributions that make up the baseline measure for the dirichlet process
        std::map<baseline_tag_t,component_model_t*> baseline_models;

        // assignments to clusters
        data_t<cluster_tag_t>* _cluster_assignments;
};

std::ostream& operator<< (std::ostream& o, const mixture_state_t& state);

#endif /* MIXTURE_STATE_HH */
