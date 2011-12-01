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

#include <algorithm>
#include <sstream>
#include <string.h>
#include <math.h>

#include <component-model.hh>

using namespace std;

// Constants
////////////////////////////////////////////////////////////////////////////////

static const size_t ALPHABET_SIZE = 4; 
static const size_t CONTEXT       = 2; 

// Standard output
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const nucleotide_sequence_t& sequence)
{
        for (size_t i = 0; i < sequence.size(); i++) {
                o << sequence[i];
        }

        return o;
}

// Test functions
////////////////////////////////////////////////////////////////////////////////

static void sanity_check() {
        size_t k = ALPHABET_SIZE*ALPHABET_SIZE;
        char n[4] = {'A', 'C', 'T', 'G'};
        string n_str("AGGGGACGTCGATGCGTGATCGACTACGGCXX");
        vector<size_t> sizes(k, n_str.length());
        sequence_data_t<short> data;
        sequence_data_t<cluster_tag_t> cluster_assignments(sizes, 0);

        // generate sequences
        for (short i = 0; i < (short)ALPHABET_SIZE; i++) {
                n_str[n_str.length()-2] = n[i];
                for (short j = 0; j < (short)ALPHABET_SIZE; j++) {
                        n_str[n_str.length()-1] = n[j];
                        nucleotide_sequence_t sequence(n_str);
                        data.push_back(sequence);
                }
        }
        //cout << data << endl;

        double sum = 0;
        for (size_t i = 0; i < k; i++) {
                cluster_assignments[i][n_str.length()-1] = 2;
                cluster_assignments[i][n_str.length()-2] = 2;
        }
        //cout << cluster_assignments << endl;
        MarkovChainMixture model(ALPHABET_SIZE, CONTEXT, data, cluster_assignments, 0);
        for (size_t i = 0; i < k; i++) {
                cluster_assignments[i][n_str.length()-1] = 1;
                cluster_assignments[i][n_str.length()-2] = 1;
        }
        model.add(range_t(seq_index_t(0,0), n_str.length()-2));
        for (size_t i = 0; i < k; i++) {
                double p = model.predictive(range_t(seq_index_t(i,n_str.length()-2), 2));
                sum += p;
                cout << *static_cast<const nucleotide_sequence_t*>(&data[i])
                     << ": "
                     << p
                     << endl;
        }

        // fill counts1 with sequences
        cout << "Sum: " << sum << endl;
}

int main()
{
        sanity_check();

        return 0;
}
