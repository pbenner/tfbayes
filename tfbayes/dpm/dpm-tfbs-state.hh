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
         dpm_tfbs_state_t(const tfbs_options_t& options,
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
        dpm_tfbs_state_t& operator=(const dpm_tfbs_state_t& state);

        // methods
        ////////////////////////////////////////////////////////////////////////
        void add   (const range_t& range, cluster_tag_t tag);
        void remove(const range_t& range);

        bool valid_foreground_position(const range_t& range) const;
        bool get_free_range(const index_t& index, size_t& length);

        bool move_left (cluster_t& cluster, cluster_tag_t bg_cluster_tag, size_t n = 1);
        bool move_right(cluster_t& cluster, cluster_tag_t bg_cluster_tag, size_t n = 1);

        void save(cluster_tag_t cluster_tag, cluster_tag_t bg_cluster_tag);
        void restore();

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
        bool is_tfbs_start_position(const index_t& index) const;
        bool is_background(const index_t& index) const;
        bool is_background(cluster_tag_t tag) const;
        bool is_background(const cluster_t& cluster) const;
        cluster_tag_t add_background_cluster(component_model_t& component_model);
        bool set_length(cluster_t& cluster, cluster_tag_t bg_cluster_tag, size_t n);

        // data
        ////////////////////////////////////////////////////////////////////////

        // record start positions of tfbs
        sequence_data_t<short> tfbs_start_positions;

        // keep track of the number of transcription factor binding sites
        size_t num_tfbs;
        // minimum and maximum lengths of tfbs
        size_t min_foreground_length;
        size_t max_foreground_length;

        // auxiliary variables to save the current state
        dpm_tfbs_state_t* state_p;

        const data_tfbs_t* m_data;

        // the bg_cluster_tag is determined after the state is
        // initialize, so this can't be a constant
        typedef std::vector<cluster_tag_t> bg_cluster_tags_t;
        bg_cluster_tags_t bg_cluster_tags;
};

#endif /* __TFBAYES_DPM_DPM_TFBS_STATE_HH__ */
