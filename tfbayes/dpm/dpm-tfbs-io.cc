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

#include <tfbayes/dpm/dpm-tfbs-io.hh>

using namespace std;

static
ostream& operator<< (ostream& o, const matrix<double>& m)
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

void dpm_tfbs_save_result(ostream& file, dpm_tfbs_pmcmc_t& sampler)
{
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
        file << "temperature =" << endl
             << history.temperature;
        file << "partitions =" << endl
             << history.partitions;
}
