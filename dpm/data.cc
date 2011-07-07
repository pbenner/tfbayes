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

#include <clustermanager.hh>
#include <data.hh>
#include <statistics.hh>

using namespace std;

Data::Data(size_t n, char *sequences[], cluster_tag_t default_tag)
        : n_sequences(n)
{
        for(size_t i = 0; i < n; i++) {
                size_t m = strlen(sequences[i]);
                this->sequences.push_back(sequences[i]);
                this->sequences_length.push_back(m);
                this->cluster_assignments.push_back(vector<cluster_tag_t>(m, default_tag));
                for(size_t j = 0; j < m; j++) {
                        element_t e = {i, j};
                        elements.push_back(e);
                }
        }
        for (Data::iterator it = begin(); it != end(); it++) {
                elements_randomized.push_back(&(*it));
        }
        shuffle();
}

Data::~Data() {
}

void
Data::shuffle() {
        random_shuffle(elements_randomized.begin(), elements_randomized.end());
}

const element_t&
Data::operator[](size_t i) const {
        return elements[i];
}

element_t&
Data::operator[](size_t i) {
        return elements[i];
}

bool
Data::valid_for_sampling(const element_t& element, size_t length, word_t& word)
{
        const size_t sequence   = element.sequence;
        const size_t position   = element.position;
        const cluster_tag_t tag = cluster_assignments[sequence][position];

        // check if there is enough space
        if (sequences_length[sequence] - position < length) {
                return false;
        }
        // check if all successive elements belong to the same class
        for (size_t i = 1; i < length; i++) {
                if (cluster_assignments[sequence][position+i] != tag) {
                        return false;
                }
        }
        word.sequence  = sequence;
        word.position  = position;
        word.length    = length;
        word.sequences = &sequences;
        return true;
}

// ostream&
// operator<< (ostream& o, element_t const& element) {
//         o << element.tag << ":(" << element.x[0] << ","
//           << element.x[1] << ")@" << element.original_cluster;

//         return o;
// }

// ostream&
// operator<< (ostream& o, Data const& data) {
//         for (Data::const_iterator it = data.begin(); it != data.end(); it++) {
//                 if (it != data.begin()) {
//                         o << ", ";
//                 }
//                 o << *it;
//         }

//         return o;
// }
