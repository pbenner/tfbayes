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

#ifndef TFBS_DATA_HH
#define TFBS_DATA_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>

#include "data.hh"

using namespace std;

class TfbsData : public Data {
public:
        TfbsData(int n, int m, char *sequences[], int *clusters[]);
        ~TfbsData();

        friend ostream& operator<< (std::ostream& o, TfbsData const& data);

        char get_nucleotide(const Data::element& e) const {
                return sequences[e.x[0]][e.x[1]];
        }

private:
        vector<string> sequences;
        size_t n_sequences;
        size_t sequence_length;
};

#endif /* TFBS_DATA_HH */
