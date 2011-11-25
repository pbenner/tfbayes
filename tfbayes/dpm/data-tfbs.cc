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

static inline
bool valid_sampling_index(const vector<string>& sequences, const index_t& index, size_t tfbs_length)
{
        if (sequences[index[0]].size() - index[1] < tfbs_length) {
                return false;
        }
        for (size_t i = 0; i < tfbs_length; i++) {
                if (sequences[index[0]][index[1]+i] == 'N') {
                        return false;
                }
        }
        return true;
}

DataTFBS::DataTFBS(const vector<string>& sequences, size_t tfbs_length)
        : sequence_data_t<short>(code_sequences(sequences)),
          sequences(sequences), _n_sequences(sequences.size()), _elements(0)
{
        for(size_t i = 0; i < sequences.size(); i++) {
                // store length of sequences
                this->sequences_length.push_back(sequences[i].size());
                // and a list of indices
                for(size_t j = 0; j < sequences[i].size(); j++) {
                        if (sequences[i][j] != 'N') {
                                this->indices.push_back(new seq_index_t(i,j));
                                _elements++;
                        }
                }
        }
        // generate a randomized list of indices
        for (DataTFBS::iterator it = begin(); it != end(); it++) {
                if (valid_sampling_index(sequences, **it, tfbs_length)) {
                        sampling_indices.push_back(*it);
                }
        }
        shuffle();
}

DataTFBS::DataTFBS(const DataTFBS& data)
        : sequences(data.sequences), sequences_length(data.sequences_length),
          _n_sequences(data._n_sequences), _elements(data._elements)
{
        for (DataTFBS::const_iterator it = data.begin(); it != data.end(); it++) {
                index_t* index = (**it).clone();
                indices.push_back(index);
                sampling_indices.push_back(index);
        }
        shuffle();
}

DataTFBS::~DataTFBS()
{
        for (DataTFBS::iterator it = begin(); it != end(); it++) {
                delete(*it);
        }
}

void
DataTFBS::shuffle() {
        random_shuffle(sampling_indices.begin(), sampling_indices.end());
}

const index_t&
DataTFBS::operator[](size_t i) const {
        return *indices[i];
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
