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

#ifndef __TFBAYES_DPM_DPM_TFBS_STATE_HH__
#define __TFBAYES_DPM_DPM_TFBS_STATE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <set>
#include <vector>

#include <boost/random/mersenne_twister.hpp>

#include <tfbayes/dpm/cluster.hh>
#include <tfbayes/dpm/data-tfbs.hh>
#include <tfbayes/dpm/mixture-state.hh>
#include <tfbayes/dpm/state.hh>

class dpm_tfbs_state_t : public gibbs_state_t {
public:
         dpm_tfbs_state_t(const std::vector<size_t>& sizes,
                          size_t tfbs_length,
                          const data_tfbs_t& data);
         dpm_tfbs_state_t(const dpm_tfbs_state_t& state);
        ~dpm_tfbs_state_t();

        friend void swap(dpm_tfbs_state_t& first, dpm_tfbs_state_t& second);

        virtual dpm_tfbs_state_t* clone() const;

        // auxiliary types
        ////////////////////////////////////////////////////////////////////////
        typedef cluster_t::const_iterator cl_iterator;

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual dpm_tfbs_state_t& operator=(const mixture_state_t& state);

        // methods
        ////////////////////////////////////////////////////////////////////////
        void add   (const range_t& range, cluster_tag_t tag);
        void remove(const range_t& range, cluster_tag_t tag);
        void remove(const index_i& index, cluster_tag_t tag);

        bool move_left (cluster_t& cluster, cluster_tag_t bg_cluster_tag);
        bool move_right(cluster_t& cluster, cluster_tag_t bg_cluster_tag);

        bool proposal(cluster_t& cluster, std::stringstream& ss, boost::random::mt19937& gen);
        void save(cluster_tag_t cluster_tag, cluster_tag_t bg_cluster_tag);
        void restore();

        bool valid_tfbs_position(const index_i& index) const;

        // access to cluster assignments, the data type might be
        // different in child classes, so make this virtual
        virtual sequence_data_t<cluster_tag_t>& cluster_assignments() {
                return *static_cast<sequence_data_t<cluster_tag_t>*>(&gibbs_state_t::cluster_assignments());
        }
        virtual const sequence_data_t<cluster_tag_t>& cluster_assignments() const {
                return *static_cast<const sequence_data_t<cluster_tag_t>*>(&gibbs_state_t::cluster_assignments());
        }

        dpm_partition_t partition() const;
        void set_partition(const dpm_partition_t& partition);
        bool is_background(cluster_tag_t tag) const;
        bool is_background(const cluster_t& cluster) const;
        cluster_tag_t add_background_cluster(component_model_t& component_model);

        // data
        ////////////////////////////////////////////////////////////////////////

        // record start positions of tfbs
        sequence_data_t<short> tfbs_start_positions;

        // keep track of the number of transcription factor binding sites
        size_t num_tfbs;

        // auxiliary variables to save the current state
        size_t num_tfbs_p;
        sequence_data_t<cluster_tag_t> cluster_assignments_p;
        sequence_data_t<short> tfbs_start_positions_p;
        cluster_t* cluster_p;
        cluster_t* cluster_bg_p;

        const data_tfbs_t* _data;

        // constants
        size_t tfbs_length;
        // the bg_cluster_tag is determined after the state is
        // initialize, so this can't be a constant
        typedef std::vector<cluster_tag_t> bg_cluster_tags_t;
        bg_cluster_tags_t bg_cluster_tags;
};

#endif /* __TFBAYES_DPM_DPM_TFBS_STATE_HH__ */
