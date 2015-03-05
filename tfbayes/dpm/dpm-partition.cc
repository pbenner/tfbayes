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

#include <tfbayes/dpm/dpm-partition.hh>

using namespace std;

ostream& operator<< (ostream& o, const dpm_subset_t& dpm_subset)
{
        for (dpm_subset_t::const_iterator is = dpm_subset.begin();
             is != dpm_subset.end(); is++)
        {
                if (is != dpm_subset.begin()) {
                        o << ", ";
                }
                o << is->index()
                  << ":" << is->length();
                if (is->reverse()) {
                        o << "!";
                }
        }
        return o;
}

ostream& operator<< (ostream& o, const model_id_t& id)
{
        o << id.name << ":" << id.length;
        return o;
}

ostream& operator<< (ostream& o, const dpm_partition_t& partition)
{
        for (dpm_partition_t::const_iterator it = partition.begin();
             it != partition.end(); it++) {
                if (it != partition.begin()) {
                        o << ", ";
                }
                o << it->model_id() << ":{"
                  << *it << "}";
        }
        return o;
}

ostream& operator<< (ostream& o, const dpm_partition_list_t& partitions)
{
        if (partitions.size() == 0) {
                return o;
        }
        for (size_t i = 0; i < partitions.size(); i++) {
                o << "\t" << partitions[i] << endl;
        }
        return o;
}
