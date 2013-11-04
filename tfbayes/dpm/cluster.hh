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

#ifndef __TFBAYES_DPM_CLUSTER_HH__
#define __TFBAYES_DPM_CLUSTER_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <vector>
#include <iostream>

#include <boost/unordered_set.hpp> 

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/dpm/datatypes.hh>
#include <tfbayes/dpm/observer.hh>
#include <tfbayes/utility/clonable.hh>

////////////////////////////////////////////////////////////////////////////////
// This class represents a single cluster for the dirichlet process
// mixture. Each cluster is linked to a probability distribution and
// keeps track of the sufficient statistics assiciated to each cluster.
////////////////////////////////////////////////////////////////////////////////

class cluster_t : public Observed<cluster_event_t> {
public:
         cluster_t(component_model_t* model, cluster_tag_t cluster_tag, baseline_tag_t baseline_tag,
                   bool destructible = true, bool record = false);
         cluster_t(component_model_t* model, cluster_tag_t cluster_tag, baseline_tag_t baseline_tag,
                   Observer<cluster_event_t>* observer, bool destructible = true, bool record = false);
         cluster_t(const cluster_t& cluster);
        ~cluster_t();

        // types
        typedef boost::unordered_set<range_t> elements_t;

        typedef elements_t::iterator iterator;
        typedef elements_t::const_iterator const_iterator;

        // iterators
        const_iterator begin() const { return _elements.begin(); }
        const_iterator end()   const { return _elements.end();   }

        // friends
        friend std::ostream& operator<< (std::ostream& o, const cluster_t& cluster);

        // operators
        cluster_t& operator=(const cluster_t& cluster);

        // methods
        void add_observations(const range_t& range);
        void remove_observations(const range_t& range);
        size_t size() const;
        cluster_tag_t cluster_tag() const;
        baseline_tag_t baseline_tag() const;
        bool destructible() const;
        component_model_t& model();
        const component_model_t& model() const;
        elements_t elements() const;

private:
        component_model_t* _model;
        const cluster_tag_t _cluster_tag;
        const baseline_tag_t _baseline_tag;
        const bool _destructible;
        const bool _record;
        elements_t _elements;

        // number of elements in the cluster
        size_t _size;
};

#endif /* __TFBAYES_DPM_CLUSTER_HH__ */
