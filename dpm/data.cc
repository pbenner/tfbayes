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

using namespace std;

#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include "data.hh"
#include "statistics.hh"

Data::Data(int n, int m, char *sequences[], int *clusters[])
        : n_sequences(n), sequence_length(m)
{
        int tag = 0;

        for(int i = 0; i < n; i++) {
                this->sequences.push_back(sequences[i]);
                for(int j = 0; j < m; j++) {
                        Data::x_t x;
                        x.push_back(i);
                        x.push_back(j);
                        Data::element e = {
                                x, tag++,
                                clusters[i][j]
                        };
                        elements.push_back(e);
                }
        }
        for (Data::iterator it = this->begin(); it != this->end(); it++) {
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

const Data::element&
Data::operator[](size_t i) const {
        return elements[i];
}

Data::element&
Data::operator[](size_t i) {
        return elements[i];
}

Data::iterator
Data::find(const element& elem) {
        for (iterator it = begin(); it != end(); it++) {
                if ((*it).tag == elem.tag)
                        return it;
        }
        return end();
}

char
Data::get_nucleotide(const Data::element& e) const {
        return sequences[e.x[0]][e.x[1]];
}

void
Data::get_nucleotide(const Data::element& e, int n, char *buf) const {
        for (int i = 0; i < n; i++) {
                if (e.x[1]+i < sequence_length) {
                        buf[i] = sequences[e.x[0]][e.x[1]+i];
                }
                else {
                        buf[i] = '\0';
                }
        }
}

int
Data::num_successors(const Data::element& e) {
        return sequence_length-e.x[1]-1;
}

ostream&
operator<< (ostream& o, Data::element const& element) {
        o << element.tag << ":(" << element.x[0] << ","
          << element.x[1] << ")@" << element.original_cluster;

        return o;
}

ostream&
operator<< (ostream& o, Data const& data) {
        for (Data::const_iterator it = data.begin(); it != data.end(); it++) {
                if (it != data.begin()) {
                        o << ", ";
                }
                o << *it;
        }

        return o;
}
