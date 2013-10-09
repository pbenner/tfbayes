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

#ifndef DPM_TFBS_STATE_HH
#define DPM_TFBS_STATE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <tfbayes/dpm/cluster.hh>
#include <tfbayes/dpm/data-tfbs.hh>
#include <tfbayes/dpm/mixture-state.hh>
#include <tfbayes/dpm/state.hh>

class dpm_tfbs_state_t : public gibbs_state_t {
public:
         dpm_tfbs_state_t(const std::vector<size_t>& sizes,
                          size_t tfbs_length,
                          cluster_tag_t bg_cluster_tag,
                          const data_tfbs_t& data);
         dpm_tfbs_state_t(const dpm_tfbs_state_t& state);
        ~dpm_tfbs_state_t();

        virtual dpm_tfbs_state_t* clone() const;

        // auxiliary types
        ////////////////////////////////////////////////////////////////////////
        typedef cluster_t::const_iterator cl_iterator;

        // methods
        ////////////////////////////////////////////////////////////////////////
        void add   (const index_i& index, cluster_tag_t tag);
        void remove(const index_i& index, cluster_tag_t tag);

        bool move_left (cluster_t& cluster);
        bool move_right(cluster_t& cluster);

        bool proposal(cluster_t& cluster, std::stringstream& ss);
        void save(cluster_tag_t cluster_tag);
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
        cluster_tag_t bg_cluster_tag;
};

#endif /* DPM_TFBS_STATE_HH */
