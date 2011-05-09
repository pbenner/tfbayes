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

#include "tfbs-data.hh"
#include "statistics.hh"

TfbsData::TfbsData(int n, int m, char *sequences[], int *clusters[])
        : Data(), n_sequences(n), sequence_length(m)
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
                                clusters[i][j] };
                        elements.push_back(e);
                }
        }
        for (Data::iterator it = this->begin(); it != this->end(); it++) {
                elements_randomized.push_back(&(*it));
        }
        shuffle();
}

TfbsData::~TfbsData() {
}

ostream& operator<< (ostream& o, Data::element const& element) {
        o << element.tag << ":(" << element.x[0] << ","
          << element.x[1] << ")@" << element.original_cluster;

        return o;
}

ostream& operator<< (ostream& o, TfbsData const& data) {
        for (Data::const_iterator it = data.begin(); it != data.end(); it++) {
                if (it != data.begin()) {
                        o << ", ";
                }
                o << *it;
        }

        return o;
}
