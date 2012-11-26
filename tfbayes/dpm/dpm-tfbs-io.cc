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

#include <dpm-tfbs-io.hh>

using namespace std;

void dpm_tfbs_save_result(ostream& file, dpm_tfbs_pmcmc_t& sampler)
{
        const samples_t& samples          = sampler.samples();
        const sampling_history_t& history = sampler.sampling_history();
        const double sampling_steps       = sampler.sampling_steps();

        file.setf(ios::showpoint);

        file << "[Result]" << endl;
        file << "posterior =" << endl;
        for (size_t i = 0; i < samples.probabilities.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < samples.probabilities[i].size(); j++) {
                        file << (float)samples.probabilities[i][j] << " ";
                }
                file << endl;
        }
        file << "components =" << endl;
        for (size_t i = 0; i < history.components.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < history.components[i].size(); j++) {
                        file << history.components[i][j] << " ";
                }
                file << endl;
        }
        file << "switches =" << endl;
        for (size_t i = 0; i < history.switches.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < history.switches[i].size(); j++) {
                        file << history.switches[i][j] << " ";
                }
                file << endl;
        }
        file << "likelihood =" << endl;
        for (size_t i = 0; i < history.likelihood.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < history.likelihood[i].size(); j++) {
                        file << history.likelihood[i][j] << " ";
                }
                file << endl;
        }
        file << "posterior =" << endl;
        for (size_t i = 0; i < history.posterior.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < history.posterior[i].size(); j++) {
                        file << history.posterior[i][j] << " ";
                }
                file << endl;
        }
        file << "graph = ";
        for (dpm_graph_t::const_iterator it = samples.graph.begin();
             it != samples.graph.end(); it++) {
                file << *static_cast<const seq_index_t*>(&(*it).first.index1) << "-"
                     << *static_cast<const seq_index_t*>(&(*it).first.index2) << "="
                     << static_cast<double>((*it).second)/sampling_steps << " ";
        }
        file << endl;
        file << "map_partition = ";
        for (dpm_partition_t::const_iterator it = samples.map_partition.begin();
             it != samples.map_partition.end(); it++) {
                if (it == samples.map_partition.begin()) {
                        file << "{";
                }
                else {
                        file << ", {";
                }
                for (index_set_t::const_iterator is = it->begin(); is != it->end(); is++)
                {
                        if (is != it->begin()) {
                                file << ", ";
                        }

                        file << *static_cast<const seq_index_t*>(*is);
                }
                file << "}";
        }
        file << endl;
}
