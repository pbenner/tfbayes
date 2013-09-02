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

dpm_partition_t
dpm_tfbs_median(const vector<dpm_partition_t>& partitions,
                const vector<size_t>& sizes, size_t tfbs_length)
{

        sequence_data_t<cluster_tag_t>(sizes, -1);

        return dpm_partition_t();
}
