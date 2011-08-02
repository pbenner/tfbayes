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

#include <iostream>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <ctime>

#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include <data-tfbs.hh>

using namespace std;

DataTFBS::DataTFBS(const vector<string>& sequences)
        : sequence_data_t<short>(code_sequences(sequences)), sequences(sequences), _n_sequences(sequences.size()), _elements(0)
{
        for(size_t i = 0; i < sequences.size(); i++) {
                // store length of sequences
                this->sequences_length.push_back(sequences[i].size());
                // and a list of indices
                for(size_t j = 0; j < sequences[i].size(); j++) {
                        this->indices.push_back(index_t(i,j));
                        _elements++;
                }
        }
        // generate a randomized list of indices
        for (DataTFBS::iterator it = begin(); it != end(); it++) {
                indices_randomized.push_back(&(*it));
        }
        shuffle();
}

void
DataTFBS::shuffle() {
        random_shuffle(indices_randomized.begin(), indices_randomized.end());
}

const index_t&
DataTFBS::operator[](size_t i) const {
        return indices[i];
}

size_t
DataTFBS::length() const {
        return _n_sequences;
}

size_t
DataTFBS::length(size_t i) const {
        if (i < _n_sequences) {
                return sequences_length[i];
        }
        else {
                return 0;
        }
}

const std::vector<size_t>&
DataTFBS::lengths() const
{
        return sequences_length;
}

ostream& operator<< (ostream& o, const DataTFBS& data)
{
        for (size_t i = 0; i < data.sequences.size(); i++) {
                o << data.sequences[i] << endl;
        }

        return o;
}
