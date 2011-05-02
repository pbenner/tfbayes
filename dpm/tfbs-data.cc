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
        : Data(), k(1)
{
        int tag = 0;

        for(int i = 0; i < n; i++) {
                for(int j = 0; j < m; j++) {
                        Data::element e = {
                                sequences[i][j],
                                tag++,
                                clusters[i][j] };
                        elements.push_back(e);
                }
                Data::element e = {
                        '\0', 0, 0 };
                elements.push_back(e);
        }
}

TfbsData::~TfbsData() {
}
