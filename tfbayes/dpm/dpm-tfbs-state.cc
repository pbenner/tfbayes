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

#include <algorithm>

#include <tfbayes/dpm/dpm-tfbs-state.hh>

using namespace std;

dpm_tfbs_state_t::dpm_tfbs_state_t(
        const tfbs_options_t& options,
        const data_tfbs_t& data)
        : gibbs_state_t(sequence_data_t<cluster_tag_t>(data.sizes(), -1)),
          // starting positions of tfbs
          tfbs_start_positions(data.sizes(), 0),
          // number of transcription factor binding sites
          num_tfbs(0),
          // auxiliary variables
          num_tfbs_p(0),
          cluster_assignments_p(data.sizes(), -1),
          tfbs_start_positions_p(data.sizes(), 0),
          cluster_p(NULL),
          cluster_bg_p(NULL),
          m_data(&data)
{
        assert(options.baseline_lengths.size() > 0);
        // find minimum and maximum lengths of the foreground model
        vector<size_t> baseline_lengths;
        for (size_t i = 0; i < options.baseline_lengths.size(); i++) {
                copy(options.baseline_lengths[i].begin(), options.baseline_lengths[i].end(),
                     back_inserter(baseline_lengths));
        }
        min_foreground_length = *min_element(baseline_lengths.begin(), baseline_lengths.end());
        max_foreground_length = *max_element(baseline_lengths.begin(), baseline_lengths.end());
}

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
          // length of tfbs
          min_foreground_length(state.min_foreground_length),
          max_foreground_length(state.max_foreground_length),
          // auxiliary variables
          num_tfbs_p(0),
          cluster_assignments_p(state.tfbs_start_positions.sizes(), -1),
          tfbs_start_positions_p(state.tfbs_start_positions.sizes(), 0),
          cluster_p(NULL),
          cluster_bg_p(NULL),
          m_data(state.m_data),
          bg_cluster_tags(state.bg_cluster_tags)
{ }

void swap(dpm_tfbs_state_t& first, dpm_tfbs_state_t& second)
{
        swap(static_cast<gibbs_state_t&>(first),
             static_cast<gibbs_state_t&>(second));
        swap(first.tfbs_start_positions,   second.tfbs_start_positions);
        swap(first.num_tfbs,               second.num_tfbs);
        swap(first.min_foreground_length,  second.min_foreground_length);
        swap(first.max_foreground_length,  second.max_foreground_length);
        swap(first.cluster_assignments_p,  second.cluster_assignments_p);
        swap(first.tfbs_start_positions_p, second.tfbs_start_positions_p);
        swap(first.cluster_p,              second.cluster_p);
        swap(first.cluster_bg_p,           second.cluster_bg_p);
        swap(first.m_data,                 second.m_data);
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
dpm_tfbs_state_t::valid_foreground_position(const range_t& range) const
{
        index_t current_index = range.index();

        // at the first position there either has to be background or
        // the beginning of a tfbs
        if (!is_background(current_index) && !is_tfbs_start_position(current_index)) {
                return false;
        }
        // check if there is no tfbs starting at later positions
        for (size_t i = 1; i < range.length(); i++) {
                current_index[1] = range.index()[1]+i;

                // check if index is out of range
                if (size_t(current_index[1]) >= (*m_data)[current_index[0]].size()) {
                        return false;
                }
                // check for a tfbs
                if (is_tfbs_start_position(current_index)) {
                        return false;
                }
        }
        return true;
}

bool
dpm_tfbs_state_t::get_free_range(const index_t& index, size_t& length)
{
        index_t current_index = index;

        // at the first position there either has to be background or
        // the beginning of a tfbs
        if (!is_background(index) && !is_tfbs_start_position(index)) {
                length = 0;
                return false;
        }
        // check if there is no tfbs starting at later positions
        for (size_t i = 1; i < max_foreground_length; i++) {
                current_index[1] = index[1]+i;

                // check if index is out of range
                if (size_t(current_index[1]) >= (*m_data)[current_index[0]].size()) {
                        length = i;
                        return i >= min_foreground_length;
                }
                // check for a tfbs
                if (is_tfbs_start_position(current_index)) {
                        length = i;
                        return i >= min_foreground_length;
                }
        }
        length = max_foreground_length;
        return true;
}

void
dpm_tfbs_state_t::add(const range_t& range, cluster_tag_t cluster_tag)
{
        if (is_background(cluster_tag)) {
                operator[](cluster_tag).add_observations(range);
        }
        else {
                // cluster of the foreground model starting at the
                // first position
                cluster_t& cluster = operator[](cluster_tag);
                // get the length of the foreground model
                size_t cluster_length = cluster.model().id().length;
                if (cluster_length < range.length()) {
                        // split range in two
                        range_t range1(range);
                        range_t range2(range);
                        range1.length()     = cluster_length;
                        range2.index ()[1] += cluster_length;
                        range2.length()    -= cluster_length;
                        // add observations to the foreground
                        // model
                        cluster.add_observations(range1);
                        // add remaining observations to the
                        // background model
                        add(range2, bg_cluster_tags[0]);
                }
                else {
                        cluster.add_observations(range);
                }
                tfbs_start_positions[range.index()] = range.reverse() ? -1 : 1;
                num_tfbs++;
        }
}

void
dpm_tfbs_state_t::remove(const range_t& range)
{
        const index_t& index = range.index();
        // cluster of the foreground model starting at the
        // first position
        cluster_tag_t cluster_tag = operator[](index);
        cluster_t&    cluster     = operator[](cluster_tag);
        // if there is background at the first position, we know that
        // the range is accurate for the background model
        if (is_background(index)) {
                cluster.remove_observations(range);
        }
        // otherwise, we might have to split the range into foreground
        // and background part
        else {
                assert(is_tfbs_start_position(range.index()));
                // get the length of the foreground model
                size_t cluster_length = cluster.model().id().length;
                // check if the length of the foreground model is
                // shorter than the range
                if (cluster_length < range.length()) {
                        // split range in two
                        range_t range1(range);
                        range_t range2(range);
                        range1.reverse()    = tfbs_start_positions[index] == -1;
                        range1.length()     = cluster_length;
                        range2.index ()[1] += cluster_length;
                        range2.length()    -= cluster_length;
                        // remove observations from the foreground
                        // model
                        cluster.remove_observations(range1);
                        // remove remaining observations
                        assert(is_background(range2.index()));
                        remove(range2);
                }
                else {
                        range_t range1(range);
                        range1.reverse()    = tfbs_start_positions[index] == -1;
                        cluster.remove_observations(range1);
                }
                // tfbs_start_position is required earliner, so reset
                // it here
                tfbs_start_positions[index] = 0;
                num_tfbs--;
        }
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
        cluster_p              = new cluster_t(operator[](   cluster_tag));
        cluster_bg_p           = new cluster_t(operator[](bg_cluster_tag));
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
                      range_t  new_range(*is);
                // one position to the left
                new_range.index()[1] -= n;
                remove(old_range);
                add   (old_range, bg_cluster_tag);
                if (new_range.index()[1] > 0
                    && valid_foreground_position(new_range)
                    && !is_tfbs_start_position(new_range.index())) {
                        remove(new_range);
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
                      range_t  new_range(*is);
                const size_t sequence_length = m_data->size(old_range.index()[0]);
                // one position to the left
                new_range.index()[1] += n;
                remove(old_range);
                add   (old_range, bg_cluster_tag);
                if (new_range.index()[1]+new_range.length() <= sequence_length
                    && valid_foreground_position(new_range)
                    && !is_tfbs_start_position(new_range.index())) {
                        remove(new_range);
                        add   (new_range, cluster.cluster_tag());
                }
        }
        return true;
}

bool
dpm_tfbs_state_t::set_length(cluster_t& cluster, cluster_tag_t bg_cluster_tag, size_t n)
{
        const cluster_t::elements_t old_elements(cluster.elements());
              cluster_t::elements_t new_elements;

        // remove elements from the cluster
        for (cl_iterator is = old_elements.begin(); is != old_elements.end(); is++) {
                const range_t range(is->index(), n, is->reverse());
                // since only the length was changed, it is not
                // necessary to check for the beginning of a tfbs
                if (range.index()[1]+range.length() <= m_data->size(range.index()[0])
                    && valid_foreground_position(range)) {
                        new_elements.insert(range);
                }
        }
        if (new_elements.begin() == new_elements.end()) {
                return false;
        }
        // remove elements
        for (cl_iterator is = old_elements.begin(); is != old_elements.end(); is++) {
                remove(*is);
                add   (*is, bg_cluster_tag);
        }
        // set the model length
        static_cast<product_dirichlet_t&>(cluster.model()).set_length(n);
        // add elements to the cluster
        for (cl_iterator is = new_elements.begin(); is != new_elements.end(); is++) {
                remove(*is);
                add   (*is, cluster.cluster_tag());
        }
        return true;
}

dpm_partition_t
dpm_tfbs_state_t::partition() const
{
        dpm_partition_t dpm_partition;

        // loop through all clusters
        for (dpm_tfbs_state_t::const_iterator it = begin(); it != end(); it++) {
                const cluster_t& cluster = **it;
                if (!is_background(cluster)) {
                        model_id_t id = cluster.model().id();
                        dpm_partition.add_component(id);
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
                                remove(*is);
                                add   (*is, bg_cluster_tags[0]);
                        }
                }
        }
        ////////////////////////////////////////////////////////////////////////////////
        // set state to given partition
        for (dpm_partition_t::const_iterator it = partition.begin(); it != partition.end(); it++) {
                const dpm_subset_t& subset(*it);
                cluster_t& cluster = get_free_cluster(subset.model_id());

                for (dpm_subset_t::const_iterator is = subset.begin(); is != subset.end(); is++) {
                        assert(valid_foreground_position(*is));
                        assert(operator[](is->index()) == bg_cluster_tags[0]);
                        remove(*is);
                        add   (*is, cluster.cluster_tag());
                }
        }
}

bool
dpm_tfbs_state_t::is_tfbs_start_position(const index_t& index) const
{
        return tfbs_start_positions[index] != 0;
}

bool
dpm_tfbs_state_t::is_background(const index_t& index) const
{
        return cluster_assignments()[index] == 0;
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

