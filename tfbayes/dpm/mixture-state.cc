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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <list>
#include <map>
#include <set>
#include <assert.h>

#include <boost/foreach.hpp>

#include <tfbayes/dpm/mixture-state.hh>

using namespace std;

mixture_state_t::mixture_state_t(const data_i<cluster_tag_t>& cluster_assignments)
        : used_clusters_size    (0)
        , free_clusters_size    (0)
        , m_cluster_assignments (cluster_assignments.clone())
{
}

mixture_state_t::~mixture_state_t() {
        // free clusters
        for (mixture_state_t::iterator_all it = clusters.begin();
             it != clusters.end(); it++) {
                delete(*it);
        }
        // free baseline models
        for (map<baseline_tag_t, component_model_t*>::iterator it = baseline_models.begin();
             it != baseline_models.end(); it++) {
                delete(it->second);
        }
        delete(m_cluster_assignments);
}

mixture_state_t::mixture_state_t(const mixture_state_t& cm)
        : used_clusters_size    (cm.used_clusters_size)
        , free_clusters_size    (cm.free_clusters_size)
        , m_cluster_assignments (cm.m_cluster_assignments->clone())
{
        for (size_t i = 0; i < cm.clusters.size(); i++) {
                cluster_t* c = new cluster_t(*cm.clusters[i]);
                c->set_observer(this);
                c->model().set_cluster_assignments(*m_cluster_assignments);
                clusters.push_back(c);
                if (c->size() == 0 && c->destructible()) {
                        free_clusters.push_back(c);
                }
                else {
                        used_clusters.push_back(c);
                }
        }
        for (map<baseline_tag_t, component_model_t*>::const_iterator it = cm.baseline_models.begin();
             it != cm.baseline_models.end(); it++) {
                baseline_models[it->first] = it->second->clone();
                baseline_models[it->first]->set_cluster_assignments(*m_cluster_assignments);
        }
}

mixture_state_t*
mixture_state_t::clone() const
{
        return new mixture_state_t(*this);
}

void swap(mixture_state_t& first, mixture_state_t& second)
{
        swap(first.clusters,              second.clusters);
        swap(first.used_clusters,         second.used_clusters);
        swap(first.free_clusters,         second.free_clusters);
        swap(first.used_clusters_size,    second.used_clusters_size);
        swap(first.free_clusters_size,    second.free_clusters_size);
        swap(first.baseline_models,       second.baseline_models);
        swap(first.m_cluster_assignments, second.m_cluster_assignments);
}

mixture_state_t&
mixture_state_t::operator=(const mixture_state_t& mixture_state)
{
        mixture_state_t tmp(mixture_state);
        swap(*this, tmp);
        return *this;
}

// add a cluster with specific model which will
// not be destructible (for other components of the mixture model)
// and not included in the baseline measures
cluster_tag_t
mixture_state_t::add_cluster(component_model_t* model)
{
        cluster_tag_t cluster_tag = clusters.size();

        // since this cluster is a fixed component of the
        // process, it should never be put into the list
        // of free clusters, so we do not need to observe it
        // and we can simply push it to the list of used
        // clusters
        cluster_t* c = new cluster_t(model, cluster_tag, -1, this, false, false);
        clusters.push_back(c);
        used_clusters.push_back(clusters.back());
        used_clusters_size++;

        return cluster_tag;
}

// add cluster with default model
cluster_tag_t
mixture_state_t::add_cluster(baseline_tag_t baseline_tag)
{
        // generate a new cluster tag
        cluster_tag_t cluster_tag = clusters.size();

        component_model_t* model  = baseline_models[baseline_tag]->clone();

        cluster_t* c = new cluster_t(model, cluster_tag, baseline_tag, this, true, true);
        clusters.push_back(c);
        // this cluster is empty, so place it into the list
        // of free clusters
        free_clusters.push_back(clusters.back());
        free_clusters_size++;

        return cluster_tag;
}

baseline_tag_t
mixture_state_t::add_baseline_model(component_model_t* distribution)
{
        baseline_tag_t baseline_tag = baseline_models.size();
        // add a baseline model to the list of distributions
        baseline_models[baseline_tag] = distribution;
        // return the new tag
        return baseline_tag;
}

cluster_t&
mixture_state_t::get_free_cluster(baseline_tag_t baseline_tag) {
        // check free clusters
        for (std::list<cluster_t*>::iterator it = free_clusters.begin();
             it != free_clusters.end(); it++) {
                if ((*it)->baseline_tag() == baseline_tag) {
                        return **it;
                }
        }
        // create new cluster
        cluster_tag_t cluster_tag = add_cluster(baseline_tag);

        return operator[](cluster_tag);
}

cluster_t&
mixture_state_t::get_free_cluster(const model_id_t& model_id) {
        // check free clusters
        for (std::list<cluster_t*>::iterator it = free_clusters.begin();
             it != free_clusters.end(); it++) {
                if ((*it)->model().id() == model_id) {
                        return **it;
                }
        }
        assert(false);
}

void
mixture_state_t::update(Observed<cluster_event_t>* observed, cluster_event_t event)
{
        cluster_t* cluster = (cluster_t*)observed;

        switch (event) {
        case cluster_event_empty:
                if (cluster->destructible()) {
                        used_clusters.remove(cluster);
                        free_clusters.push_back(cluster);
                        used_clusters_size--;
                        free_clusters_size++;
                }
                break;
        case cluster_event_nonempty:
                if (cluster->destructible()) {
                        used_clusters.push_back(cluster);
                        free_clusters.remove(cluster);
                        used_clusters_size++;
                        free_clusters_size--;
                }
                break;
        default:
                break;
        }
}

void
mixture_state_t::update(Observed<cluster_event_t>* observed, cluster_event_t event, const range_t& range)
{
        cluster_t* cluster = (cluster_t*)observed;
        iterator_t<cluster_tag_t> iterator = cluster_assignments()[range];

        switch (event) {
        case cluster_event_add_word:
                do {
                        *iterator = cluster->cluster_tag();
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

dpm_partition_t
mixture_state_t::partition() const
{
        dpm_partition_t dpm_partition;

        // loop through all clusters
        for (const_iterator it = begin(); it != end(); it++) {
                const cluster_t& cluster = **it;
                model_id_t id = cluster.model().id();
                dpm_partition.add_component(id);
                // loop through cluster elements
                for (cluster_t::const_iterator is = cluster.begin();
                     is != cluster.end(); is++) {
                        dpm_partition.back().insert(range_t(is->index(), 1));
                }
        }
        return dpm_partition;
}

ostream& operator<< (ostream& o, const mixture_state_t& state)
{
        typedef map<size_t, multiset<size_t, greater<size_t> > > map1_t;
        typedef map<std::string, map1_t>  map2_t;

        map2_t map;

        o << " -> number of clusters: " << state.size() << endl;

        BOOST_FOREACH(const cluster_t* cluster, state) {
                model_id_t id = cluster->model().id();

                map[id.name][id.length].insert(
                        cluster->size());
        }
        for (map2_t::iterator it = map.begin(); it != map.end(); it++) {

                o << boost::format(" clusters with model `%s\':") % it->first
                  << endl;
                for (map1_t::iterator is = it->second.begin(); is != it->second.end(); is++) {
                        o << boost::format(" -> and length %d have number of elements: ") % is->first;
                        BOOST_FOREACH(const size_t& length, is->second) {
                                o << length << " ";
                        }
                        o << endl;
                }
        }
        return o;
}
