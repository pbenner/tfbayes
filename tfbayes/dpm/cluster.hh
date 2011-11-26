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

#ifndef CLUSTER_HH
#define CLUSTER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <vector>
#include <iostream>

#include <clonable.hh>
#include <component-model.hh>
#include <datatypes.hh>
#include <observer.hh>

////////////////////////////////////////////////////////////////////////////////
// This class represents a single cluster for the dirichlet process
// mixture. Each cluster is linked to a probability distribution and
// keeps track of the sufficient statistics assiciated to each cluster.
////////////////////////////////////////////////////////////////////////////////

class Cluster : public Observed<cluster_event_t> {
public:
         Cluster(ComponentModel* model, cluster_tag_t cluster_tag, model_tag_t model_tag);
         Cluster(ComponentModel* model, cluster_tag_t cluster_tag, model_tag_t model_tag, bool destructible);
         Cluster(ComponentModel* model, cluster_tag_t cluster_tag, model_tag_t model_tag, Observer<cluster_event_t>* observer);
         Cluster(ComponentModel* model, cluster_tag_t cluster_tag, model_tag_t model_tag, Observer<cluster_event_t>* observer, bool destructible);
         Cluster(const Cluster& cluster);
        ~Cluster();

        // friends
        friend std::ostream& operator<< (std::ostream& o, const Cluster& cluster);

        // methods
        void add_observations(const range_t& range);
        void remove_observations(const range_t& range);
        size_t size() const;
        cluster_tag_t cluster_tag() const;
        model_tag_t model_tag() const;
        bool destructible() const;
        ComponentModel& model();

private:
        ComponentModel* _model;
        cluster_tag_t _cluster_tag;
        model_tag_t _model_tag;
        bool _destructible;
        // number of elements in the cluster
        size_t _size;
};

#endif /* CLUSTER_HH */
