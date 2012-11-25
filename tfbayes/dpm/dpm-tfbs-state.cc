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

#include <sstream>

#include <dpm-tfbs-state.hh>

using namespace std;

dpm_tfbs_state_t::dpm_tfbs_state_t(
        const std::vector<size_t>& sizes,
        size_t tfbs_length,
        cluster_tag_t bg_cluster_tag,
        const data_tfbs_t& data)
        : state_t(cluster_assignments),
          hybrid_state_t(cluster_assignments),
          cluster_assignments(sizes, -1),
          // starting positions of tfbs
          tfbs_start_positions(sizes, 0),
          // number of transcription factor binding sites
          num_tfbs(0),
          // auxiliary variables
          num_tfbs_p(0),
          cluster_assignments_p(sizes, -1),
          tfbs_start_positions_p(sizes, 0),
          cluster_p(NULL),
          cluster_bg_p(NULL),
          _data(data),
          tfbs_length(tfbs_length),
          bg_cluster_tag(bg_cluster_tag)
{ }

dpm_tfbs_state_t::~dpm_tfbs_state_t() {
        if (cluster_p != NULL) {
                delete(cluster_p);
        }
        if (cluster_bg_p != NULL) {
                delete(cluster_bg_p);
        }
}

dpm_tfbs_state_t*
dpm_tfbs_state_t::clone() const {
        return new dpm_tfbs_state_t(*this);
}

bool
dpm_tfbs_state_t::valid_tfbs_position(const index_i& index) const
{
        const size_t sequence = index[0];
        const size_t position = index[1];

        // check if there is a tfbs starting here, if not check
        // succeeding positions
        if (tfbs_start_positions[index] == 0) {
                // check if this element belongs to a tfbs that starts
                // earlier in the sequence
                if (operator[](index) != bg_cluster_tag) {
                        return false;
                }
                if (cluster_assignments[seq_index_t(sequence, position)] == -1) {
                        return false;
                }
                // check if there is a tfbs starting within the word
                for (size_t i = 1; i < tfbs_length; i++) {
                        if (tfbs_start_positions[seq_index_t(sequence, position+i)] == 1) {
                                return false;
                        }
                        if (cluster_assignments[seq_index_t(sequence, position+i)] == -1) {
                                return false;
                        }
                }
        }

        return true;
}

void
dpm_tfbs_state_t::add(const index_i& index, cluster_tag_t tag)
{
        const range_t range(index, tfbs_length);

        if (tag == bg_cluster_tag) {
                operator[](tag).add_observations(range);
        }
        else {
                operator[](tag).add_observations(range);
                num_tfbs++;
                tfbs_start_positions[index] = 1;
        }
}

void
dpm_tfbs_state_t::remove(const index_i& index, cluster_tag_t tag)
{
        const range_t range(index, tfbs_length);

        if (tag == bg_cluster_tag) {
                operator[](tag).remove_observations(range);
        }
        else {
                operator[](tag).remove_observations(range);
                num_tfbs--;
                tfbs_start_positions[index] = 0;
        }
}

void
dpm_tfbs_state_t::save(cluster_tag_t cluster_tag) {
        if (cluster_p != NULL) {
                delete(cluster_p);
        }
        if (cluster_bg_p != NULL) {
                delete(cluster_bg_p);
        }
        num_tfbs_p = num_tfbs;
        cluster_assignments_p  = cluster_assignments;
        tfbs_start_positions_p = tfbs_start_positions;
        cluster_p    = new cluster_t(operator[](cluster_tag));
        cluster_bg_p = new cluster_t(operator[](bg_cluster_tag));
}

void
dpm_tfbs_state_t::restore() {
        num_tfbs = num_tfbs_p;
        cluster_assignments  = cluster_assignments_p;
        tfbs_start_positions = tfbs_start_positions_p;
        operator[](cluster_p->cluster_tag())    = *cluster_p;
        operator[](cluster_bg_p->cluster_tag()) = *cluster_bg_p;
}

bool
dpm_tfbs_state_t::move_left(cluster_t& cluster)
{
        const cluster_t::elements_t elements(cluster.elements());

        for (cl_iterator is = elements.begin(); is != elements.end(); is++) {
                if (cluster.size() == 1) {
                        restore();
                        return false;
                }
                const range_t& range  = *is;
                const size_t sequence = range.index[0];
                const size_t position = range.index[1];
                const seq_index_t index(sequence, position-1);
                remove(range.index, cluster.cluster_tag());
                add(range.index, bg_cluster_tag);
                if (position > 0 && valid_tfbs_position(index)) {
                        remove(index, bg_cluster_tag);
                        add(index, cluster.cluster_tag());
                }
        }

        return true;
}

bool
dpm_tfbs_state_t::move_right(cluster_t& cluster)
{
        const cluster_t::elements_t elements(cluster.elements());

        for (cl_iterator is = elements.begin(); is != elements.end(); is++) {
                if (cluster.size() == 1) {
                        restore();
                        return false;
                }
                const range_t range(*is);
                const size_t sequence = range.index[0];
                const size_t position = range.index[1];
                const size_t sequence_length = _data.size(sequence);
                const seq_index_t index(sequence, position+1);
                remove(range.index, cluster.cluster_tag());
                add(range.index, bg_cluster_tag);
                if (position+tfbs_length < sequence_length && valid_tfbs_position(index)) {
                        remove(index, bg_cluster_tag);
                        add(index, cluster.cluster_tag());
                }
        }

        return true;
}

bool
dpm_tfbs_state_t::proposal(cluster_t& cluster, stringstream& ss)
{
        save(cluster.cluster_tag());

        if (cluster.cluster_tag() != bg_cluster_tag && cluster.size() > 1) {
                if (rand() % 2 == 0) {
                        ss << "move to right";
                        return move_right(cluster);
                }
                else {
                        ss << "move to left";
                        return move_left(cluster);
                }
        }
        return false;
}

#include <graph.hh>

mixture_partition_t
dpm_tfbs_state_t::mixture_partition() const
{
        mixture_partition_t mixture_partition;

        // loop through all clusters
        for (const_iterator it = this->begin(); it != this->end(); it++) {
                const cluster_t& cluster = **it;
                if (cluster.cluster_tag() != bg_cluster_tag) {
                        mixture_partition.add_component();
                        // loop through cluster elements
                        for (cl_iterator is = cluster.begin(); is != cluster.end(); is++) {
                                mixture_partition.back().insert(is->index);
                        }
                }
        }
        return mixture_partition;
}
