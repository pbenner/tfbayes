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

#ifndef STATE_HH
#define STATE_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <ostream>

#include <mixture-state.hh>

class state_t : public mixture_state_t {
public:
        state_t(data_t<cluster_tag_t>& cluster_assignments)
                : mixture_state_t(cluster_assignments) {}

        virtual ~state_t() {}

        virtual void print(std::ostream& o) const = 0;
};

class gibbs_state_t : virtual public state_t {
public:
        gibbs_state_t(data_t<cluster_tag_t>& cluster_assignments)
                : state_t(cluster_assignments) {}

        virtual void add(const index_i& index, cluster_tag_t tag) = 0;
        virtual void remove(const index_i& index, cluster_tag_t tag) = 0;
};

class metropolis_state_t : virtual public state_t {
public:
        metropolis_state_t(data_t<cluster_tag_t>& cluster_assignments)
                : state_t(cluster_assignments) {}

        virtual bool proposal(Cluster& cluster) = 0;
        virtual void restore() = 0;
};

class hybrid_state_t : public gibbs_state_t, public metropolis_state_t {
public:
        hybrid_state_t(data_t<cluster_tag_t>& cluster_assignments)
                : state_t(cluster_assignments),
                  gibbs_state_t(cluster_assignments),
                  metropolis_state_t(cluster_assignments)
                {}
};

std::ostream& operator<< (std::ostream& o, const gibbs_state_t& state);

#endif /* STATE_HH */
