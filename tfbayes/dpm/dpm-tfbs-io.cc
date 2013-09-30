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

#include <tfbayes/dpm/dpm-tfbs-io.hh>

using namespace std;

static
ostream& operator<< (ostream& o, const vector<vector<double> >& m)
{
        for (size_t i = 0; i < m.size(); i++) {
                o << "\t";
                for (size_t j = 0; j < m[i].size(); j++) {
                        o << m[i][j] << " ";
                }
                o << endl;
        }
        return o;
}

static
ostream& operator<< (ostream& o, const vector<vector<size_t> >& m)
{
        for (size_t i = 0; i < m.size(); i++) {
                o << "\t";
                for (size_t j = 0; j < m[i].size(); j++) {
                        o << m[i][j] << " ";
                }
                o << endl;
        }
        return o;
}

static
ostream& operator<< (ostream& o, const dpm_subset_t& dpm_subset)
{
        for (dpm_subset_t::const_iterator is = dpm_subset.begin();
             is != dpm_subset.end(); is++)
        {
                if (is != dpm_subset.begin()) {
                        o << ", ";
                }
                o << *static_cast<const seq_index_t*>(*is);
        }
        return o;
}

static
ostream& operator<< (ostream& o, const dpm_partition_t& partition)
{
        for (dpm_partition_t::const_iterator it = partition.begin();
             it != partition.end(); it++) {
                if (it != partition.begin()) {
                        o << ", ";
                }
                o << it->dpm_subset_tag() << ":{"
                  << *it << "}";
        }
        return o;
}

static
ostream& operator<< (ostream& o, const std::vector<dpm_partition_t>& partitions)
{
        for (std::vector<dpm_partition_t>::const_iterator it = partitions.begin();
             it != partitions.end(); it++) {
                o << "\t" << *it << endl;
        }

        return o;
}

void dpm_tfbs_save_result(ostream& file, dpm_tfbs_pmcmc_t& sampler)
{
        const samples_t& samples          = sampler.samples();
        const sampling_history_t& history = sampler.sampling_history();

        file.setf(ios::showpoint);

        file << "[Result]" << endl;
        file << "components =" << endl
             << history.components;
        file << "switches =" << endl
             << history.switches;
        file << "likelihood =" << endl
             << history.likelihood;
        file << "posterior =" << endl
             << history.posterior;
        file << "partitions =" << endl
             << samples.partitions;
        file << "map_partition = "
             << samples.map_partition
             << endl;
}
