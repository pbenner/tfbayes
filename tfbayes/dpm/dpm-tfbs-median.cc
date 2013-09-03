/* Copyright (C) 2013 Philipp Benner
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

#include <list>

#include <tfbayes/dpm/dpm-tfbs-median.hh>
#include <tfbayes/dpm/data.hh>

using namespace std;

static void
init_data(const dpm_partition_t& partition, sequence_data_t<cluster_tag_t>& data,
          size_t tfbs_length)
{
        cluster_tag_t k = 1;

        for (dpm_partition_t::const_iterator it = partition.begin();
             it != partition.end(); it++, k++) {
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const seq_index_t& index = *static_cast<const seq_index_t*>(*is);
                        for (size_t i = 0; i < tfbs_length; i++) {
                                data[index[0]][index[1]+i] = k;
                        }
                }
        }
}

static size_t
distance(const list<seq_index_t> indices,
         const sequence_data_t<cluster_tag_t>& a,
         const sequence_data_t<cluster_tag_t>& b,
         sequence_data_t<short>& checked)
{
        size_t d = 0;

        for (list<seq_index_t>::const_iterator it = indices.begin();
             it != indices.end(); it++) {
                // make sure we don't count positions twice
                if (checked[*it]) continue;
                for (list<seq_index_t>::const_iterator is = it;
                     is != indices.end(); is++) {
                        if ((a[*it] == a[*is] && b[*it] != b[*is]) ||
                            (a[*it] != a[*is] && b[*it] == b[*is])) {
                                d += 1;
                        }
                }
                checked[*it] = true;
        }
        return d;
}

static size_t
distance(const dpm_partition_t& pi_a, const dpm_partition_t& pi_b,
         const vector<size_t>& sizes, size_t tfbs_length)
{
        size_t d = 0;

        sequence_data_t<cluster_tag_t> a(sizes, 0);
        sequence_data_t<cluster_tag_t> b(sizes, 0);
        sequence_data_t<short> checked(sizes, 0);

        init_data(pi_a, a, tfbs_length);
        init_data(pi_b, b, tfbs_length);

        // go through subsets of partition a
        for (dpm_partition_t::const_iterator it = pi_a.begin();
             it != pi_a.end(); it++) {
                list<seq_index_t> indices;
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const seq_index_t& tmp = *static_cast<const seq_index_t*>(*is);
                        for (size_t i = 0; i < tfbs_length; i++) {
                                indices.push_back(seq_index_t(tmp[0], tmp[1]+i));
                        }
                }
                d += distance(indices, a, b, checked);
        }
        // go through subsets of partition b
        for (dpm_partition_t::const_iterator it = pi_b.begin();
             it != pi_b.end(); it++) {
                list<seq_index_t> indices;
                for (dpm_subset_t::const_iterator is = it->begin();
                     is != it->end(); is++) {
                        const seq_index_t& tmp = *static_cast<const seq_index_t*>(*is);
                        for (size_t i = 0; i < tfbs_length; i++) {
                                indices.push_back(seq_index_t(tmp[0], tmp[1]+i));
                        }
                }
                d += distance(indices, a, b, checked);
        }
        return d;
}

const dpm_partition_t&
dpm_tfbs_median(const vector<dpm_partition_t>& partitions,
                const vector<size_t>& sizes, size_t tfbs_length)
{


        return partitions[0];
}
