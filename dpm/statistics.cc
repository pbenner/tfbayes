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

#include <math.h>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include "data.hh"
#include "statistics.hh"

using namespace std;

gsl_rng* _r;

ProductDirichlet::ProductDirichlet()
        : counts(NULL)
{
        update(NULL);
}

ProductDirichlet::ProductDirichlet(gsl_matrix* counts)
        : counts(NULL)
{
        update(counts);
}

ProductDirichlet::~ProductDirichlet() {
        if (this->counts) {
                gsl_matrix_free(this->counts);
        }
}

void ProductDirichlet::update(
        gsl_matrix* counts)
{
        if (this->counts) {
                gsl_matrix_free(this->counts);
        }

        this->counts = counts;
}

void
ProductDirichlet::remove_from_count_statistic(const char *nucleotides) {
        int len = counts->size1;

        for (int i = 0; i < len; i++) {
                switch (nucleotides[i]) {
                case 'A':
                case 'a':
                        gsl_matrix_set(counts, i, 0,
                                       gsl_matrix_get(counts, i, 0)-1);
                        break;
                case 'C':
                case 'c':
                        gsl_matrix_set(counts, i, 1,
                                       gsl_matrix_get(counts, i, 1)-1);
                        break;
                case 'G':
                case 'g':
                        gsl_matrix_set(counts, i, 2,
                                       gsl_matrix_get(counts, i, 2)-1);
                        break;
                case 'T':
                case 't':
                        gsl_matrix_set(counts, i, 3,
                                       gsl_matrix_get(counts, i, 3)-1);
                        break;
                }
        }
}

void
ProductDirichlet::add_to_count_statistic(const char *nucleotides) {
        int len = counts->size1;

        for (int i = 0; i < len; i++) {
                switch (nucleotides[i]) {
                case 'A':
                case 'a':
                        gsl_matrix_set(counts, i, 0,
                                       gsl_matrix_get(counts, i, 0)+1);
                        break;
                case 'C':
                case 'c':
                        gsl_matrix_set(counts, i, 1,
                                       gsl_matrix_get(counts, i, 1)+1);
                        break;
                case 'G':
                case 'g':
                        gsl_matrix_set(counts, i, 2,
                                       gsl_matrix_get(counts, i, 2)+1);
                        break;
                case 'T':
                case 't':
                        gsl_matrix_set(counts, i, 3,
                                       gsl_matrix_get(counts, i, 3)+1);
                        break;
                }
        }
}

double ProductDirichlet::pdf(char *buf) {
        double result = 1;

        for (unsigned int i = 0; i < this->counts->size1; i++) {
                double sum = 0;
                for (unsigned int j = 0; j < this->counts->size2; j++) {
                        sum += gsl_matrix_get(this->counts, i, j);
                }
                switch (buf[i]) {
                case 'A':
                case 'a':
                        result *= gsl_matrix_get(this->counts, i, 0)/sum;
                break;
                case 'C':
                case 'c':
                        result *= gsl_matrix_get(this->counts, i, 1)/sum;
                break;
                case 'G':
                case 'g':
                        result *= gsl_matrix_get(this->counts, i, 2)/sum;
                break;
                case 'T':
                case 't':
                        result *= gsl_matrix_get(this->counts, i, 3)/sum;
                break;
                }
        }
        return result;
}
