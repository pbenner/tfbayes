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
#include <assert.h>

#include <mixture-state.hh>

using namespace std;

mixture_state_t::mixture_state_t(data_t<cluster_tag_t>& cluster_assignments)
        : used_clusters_size(0), free_clusters_size(0),
          cluster_assignments(cluster_assignments)
{
}

mixture_state_t::~mixture_state_t() {
        // free clusters
        for (mixture_state_t::iterator_all it = clusters.begin();
             it != clusters.end(); it++) {
                delete(*it);
        }
        // free baseline models
        for (vector<component_model_t*>::iterator it = baseline_models.begin();
             it != baseline_models.end(); it++) {
                delete(*it);
        }
}

mixture_state_t::mixture_state_t(const mixture_state_t& cm)
        : used_clusters_size(cm.used_clusters_size),
          free_clusters_size(cm.free_clusters_size),
          cluster_assignments(cm.cluster_assignments)
{
        for (size_t i = 0; i < cm.clusters.size(); i++) {
                cluster_t* c = new cluster_t(*cm.clusters[i]);
                c->set_observer(this);
                clusters.push_back(c);
                if (c->size() == 0 && c->destructible()) {
                        free_clusters.push_back(c);
                }
                else {
                        used_clusters.push_back(c);
                }
        }
        for (vector<component_model_t*>::const_iterator it = cm.baseline_models.begin();
             it != cm.baseline_models.end(); it++) {
                baseline_models.push_back((*it)->clone());
        }
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
mixture_state_t::add_cluster(model_tag_t model_tag)
{
        assert(model_tag < (model_tag_t)baseline_models.size());

        cluster_tag_t cluster_tag = clusters.size();

        cluster_t* c = new cluster_t(baseline_models[model_tag]->clone(), cluster_tag, model_tag, this, true, true);
        clusters.push_back(c);
        // this cluster is empty, so place it into the list
        // of free clusters
        free_clusters.push_back(clusters.back());
        free_clusters_size++;

        return cluster_tag;
}

model_tag_t
mixture_state_t::add_baseline_model(component_model_t* distribution)
{
        // add a baseline model to the list of distributions
        baseline_models.push_back(distribution);

        // return its indes as model_tag
        return baseline_models.size()-1;
}

cluster_t&
mixture_state_t::get_free_cluster(model_tag_t model_tag) {
        // check free clusters
        for (std::list<cluster_t*>::iterator it = free_clusters.begin();
             it != free_clusters.end(); it++) {
                if ((*it)->model_tag() == model_tag) {
                        return **it;
                }
        }
        // create new cluster
        cluster_tag_t cluster_tag = add_cluster(model_tag);

        return (*this)[cluster_tag];
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
mixture_state_t::update(Observed<cluster_event_t>* observed, cluster_event_t event, const range_t& range)
{
        cluster_t* cluster = (cluster_t*)observed;
        iterator_t<cluster_tag_t> iterator = cluster_assignments[range];

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
