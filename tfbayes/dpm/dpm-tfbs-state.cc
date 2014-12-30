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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <sstream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <tfbayes/dpm/dpm-tfbs-state.hh>

using namespace std;

dpm_tfbs_state_t::dpm_tfbs_state_t(
        const std::vector<size_t>& sizes,
        size_t tfbs_length,
        const data_tfbs_t& data)
        : gibbs_state_t(sequence_data_t<cluster_tag_t>(sizes, -1)),
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
          _data(&data),
          tfbs_length(tfbs_length)
{ }

dpm_tfbs_state_t::~dpm_tfbs_state_t() {
        if (cluster_p != NULL) {
                delete(cluster_p);
        }
        if (cluster_bg_p != NULL) {
                delete(cluster_bg_p);
        }
}

dpm_tfbs_state_t::dpm_tfbs_state_t(const dpm_tfbs_state_t& state)
        : gibbs_state_t(state),
          tfbs_start_positions(state.tfbs_start_positions),
          num_tfbs(state.num_tfbs),
          // auxiliary variables
          num_tfbs_p(0),
          cluster_assignments_p(state.tfbs_start_positions.sizes(), -1),
          tfbs_start_positions_p(state.tfbs_start_positions.sizes(), 0),
          cluster_p(NULL),
          cluster_bg_p(NULL),
          _data(state._data),
          tfbs_length(state.tfbs_length),
          bg_cluster_tags(state.bg_cluster_tags)
{ }

void swap(dpm_tfbs_state_t& first, dpm_tfbs_state_t& second)
{
        swap(static_cast<gibbs_state_t&>(first),
             static_cast<gibbs_state_t&>(second));
        swap(first.tfbs_start_positions,   second.tfbs_start_positions);
        swap(first.num_tfbs,               second.num_tfbs);
        swap(first.cluster_assignments_p,  second.cluster_assignments_p);
        swap(first.tfbs_start_positions_p, second.tfbs_start_positions_p);
        swap(first.cluster_p,              second.cluster_p);
        swap(first.cluster_bg_p,           second.cluster_bg_p);
        swap(first._data,                  second._data);
        swap(first.tfbs_length,            second.tfbs_length);
        swap(first.bg_cluster_tags,        second.bg_cluster_tags);
}

dpm_tfbs_state_t*
dpm_tfbs_state_t::clone() const {
        return new dpm_tfbs_state_t(*this);
}

dpm_tfbs_state_t&
dpm_tfbs_state_t::operator=(const mixture_state_t& state)
{
        dpm_tfbs_state_t tmp(static_cast<const dpm_tfbs_state_t&>(state));
        swap(*this, tmp);
        return *this;
}

bool
dpm_tfbs_state_t::valid_tfbs_position(const index_i& index) const
{
        // check if there is a tfbs starting here, if not check
        // succeeding positions
        if (tfbs_start_positions[index] == 0) {
                const size_t sequence = index[0];
                const size_t position = index[1];
                // check if this element belongs to a tfbs that starts
                // earlier in the sequence
                for (size_t i = 0; i < tfbs_length; i++) {
                        if (!is_background(operator[](seq_index_t(sequence, position+i)))) {
                                return false;
                        }
                }
        }
        return true;
}

void
dpm_tfbs_state_t::add(const range_t& range, cluster_tag_t tag)
{
        operator[](tag).add_observations(range);

        if (!is_background(tag)) {
                tfbs_start_positions[range.index()] = range.reverse() ? -1 : 1;
                num_tfbs++;
        }
}

void
dpm_tfbs_state_t::remove(const range_t& range, cluster_tag_t tag)
{
        operator[](tag).remove_observations(range);

        if (!is_background(tag)) {
                assert((range.reverse() == 0 && tfbs_start_positions[range.index()] ==  1) ||
                       (range.reverse() == 1 && tfbs_start_positions[range.index()] == -1));
                num_tfbs--;
                tfbs_start_positions[range.index()] = 0;
        }
}

void
dpm_tfbs_state_t::remove(const index_i& index, cluster_tag_t tag)
{
        range_t range(index, tfbs_length, tfbs_start_positions[index] == -1);

        remove(range, tag);
}

void
dpm_tfbs_state_t::save(cluster_tag_t cluster_tag, cluster_tag_t bg_cluster_tag) {
        if (cluster_p != NULL) {
                delete(cluster_p);
        }
        if (cluster_bg_p != NULL) {
                delete(cluster_bg_p);
        }
        num_tfbs_p             = num_tfbs;
        cluster_assignments_p  = cluster_assignments();
        tfbs_start_positions_p = tfbs_start_positions;
        cluster_p    = new cluster_t(operator[](   cluster_tag));
        cluster_bg_p = new cluster_t(operator[](bg_cluster_tag));
}

void
dpm_tfbs_state_t::restore() {
        num_tfbs              = num_tfbs_p;
        cluster_assignments() = cluster_assignments_p;
        tfbs_start_positions  = tfbs_start_positions_p;
        operator[](cluster_p   ->cluster_tag()) = *cluster_p;
        operator[](cluster_bg_p->cluster_tag()) = *cluster_bg_p;
}

bool
dpm_tfbs_state_t::move_left(cluster_t& cluster, cluster_tag_t bg_cluster_tag, size_t n)
{
        const cluster_t::elements_t elements(cluster.elements());

        for (cl_iterator is = elements.begin(); is != elements.end(); is++) {
                if (cluster.size() == 1) {
                        restore();
                        return false;
                }
                const range_t& old_range(*is);
                      range_t  new_range(old_range);
                // one position to the left
                new_range.index()[1] -= n;
                remove(old_range, cluster.cluster_tag());
                add   (old_range, bg_cluster_tag);
                if (new_range.index()[1] > 0 && valid_tfbs_position(new_range.index())) {
                        remove(new_range, bg_cluster_tag);
                        add   (new_range, cluster.cluster_tag());
                }
        }

        return true;
}

bool
dpm_tfbs_state_t::move_right(cluster_t& cluster, cluster_tag_t bg_cluster_tag, size_t n)
{
        const cluster_t::elements_t elements(cluster.elements());

        for (cl_iterator is = elements.begin(); is != elements.end(); is++) {
                if (cluster.size() == 1) {
                        restore();
                        return false;
                }
                const range_t& old_range(*is);
                      range_t  new_range(old_range);
                const size_t sequence_length = _data->size(old_range.index()[0]);
                // one position to the left
                new_range.index()[1] += n;
                remove(old_range, cluster.cluster_tag());
                add   (old_range, bg_cluster_tag);
                if (new_range.index()[1]+new_range.length() <= sequence_length &&
                    valid_tfbs_position(new_range.index())) {
                        remove(new_range, bg_cluster_tag);
                        add   (new_range, cluster.cluster_tag());
                }
        }
        return true;
}

bool
dpm_tfbs_state_t::proposal(cluster_t& cluster, stringstream& ss, boost::random::mt19937& gen)
{
        if (is_background(cluster)) {
                return false;
        }
        boost::random::uniform_int_distribution<> dist   (0, 1);
        boost::random::uniform_int_distribution<> dist_bg(0, bg_cluster_tags.size()-1);
        boost::random::uniform_int_distribution<> steps  (1, cluster.elements().begin()->length()/2);
        size_t n = steps(gen);
        // select a background component at random
        cluster_tag_t bg_cluster_tag = bg_cluster_tags[dist_bg(gen)];

        save(cluster.cluster_tag(), bg_cluster_tag);

        if (dist(gen) == 0) {
                ss << boost::format("move %d steps to the right") % n;
                return move_right(cluster, bg_cluster_tag, n);
        }
        else {
                ss << boost::format("move %d steps to the left") % n;
                return move_left(cluster, bg_cluster_tag, n);
        }
}

dpm_partition_t
dpm_tfbs_state_t::partition() const
{
        dpm_partition_t dpm_partition;

        // loop through all clusters
        for (dpm_tfbs_state_t::const_iterator it = begin(); it != end(); it++) {
                const cluster_t& cluster = **it;
                if (!is_background(cluster)) {
                        dpm_partition.add_component(cluster.baseline_tag());
                        // loop through cluster elements
                        for (cl_iterator is = cluster.begin(); is != cluster.end(); is++) {
                                dpm_partition.back().insert(*is);
                        }
                }
        }
        return dpm_partition;
}

void
dpm_tfbs_state_t::set_partition(const dpm_partition_t& partition)
{
        assert(bg_cluster_tags.size() > 0);
        ////////////////////////////////////////////////////////////////////////////////
        // release all clusters
        for (dpm_tfbs_state_t::const_iterator it = begin(); it != end(); it++) {
                const cluster_t& cluster = **it;
                if (cluster.cluster_tag() != bg_cluster_tags[0]) {
                        // loop through cluster elements
                        for (cl_iterator is = cluster.begin(); is != cluster.end(); is++) {
                                remove(*is, cluster.cluster_tag());
                                add   (*is, bg_cluster_tags[0]);
                        }
                }
        }
        ////////////////////////////////////////////////////////////////////////////////
        // set state to given partition
        for (dpm_partition_t::const_iterator it = partition.begin(); it != partition.end(); it++) {
                const dpm_subset_t& subset(*it);
                cluster_t& cluster = get_free_cluster(subset.dpm_subset_tag());

                for (dpm_subset_t::const_iterator is = subset.begin(); is != subset.end(); is++) {
                        assert(valid_tfbs_position(is->index()));
                        assert(operator[](is->index()) == bg_cluster_tags[0]);
                        remove(*is, bg_cluster_tags[0]);
                        add   (*is, cluster.cluster_tag());
                }
        }
}

bool
dpm_tfbs_state_t::is_background(cluster_tag_t tag) const
{
        return find(bg_cluster_tags.begin(), bg_cluster_tags.end(), tag)
                != bg_cluster_tags.end();
}

bool
dpm_tfbs_state_t::is_background(const cluster_t& cluster) const
{
        return is_background(cluster.cluster_tag());
}

cluster_tag_t
dpm_tfbs_state_t::add_background_cluster(component_model_t& component_model)
{
        cluster_tag_t cluster_tag = add_cluster(&component_model);
        bg_cluster_tags.push_back(cluster_tag);
        return cluster_tag;
}

