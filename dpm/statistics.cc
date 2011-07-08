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

#include <data.hh>
#include <statistics.hh>

using namespace std;

gsl_rng* _r;

ProductDirichlet::ProductDirichlet(gsl_matrix* alpha)
{
        counts = gsl_matrix_alloc(alpha->size1, 4);

        for (size_t i = 0; i < alpha->size1; i++) {
                gsl_matrix_set(counts, i, 0, gsl_matrix_get(alpha, i, 0));
                gsl_matrix_set(counts, i, 1, gsl_matrix_get(alpha, i, 1));
                gsl_matrix_set(counts, i, 2, gsl_matrix_get(alpha, i, 2));
                gsl_matrix_set(counts, i, 3, gsl_matrix_get(alpha, i, 3));
        }
}

ProductDirichlet::ProductDirichlet(const ProductDirichlet& distribution)
{
        counts = gsl_matrix_alloc(distribution.counts->size1, distribution.counts->size2);
        gsl_matrix_memcpy(counts, distribution.counts);
}

ProductDirichlet::~ProductDirichlet() {
        gsl_matrix_free(counts);
}

ProductDirichlet*
ProductDirichlet::clone() const {
        return new ProductDirichlet(*this);
}

size_t
ProductDirichlet::count_observations(const word_t& word) {
        return word.length/counts->size1;
}

size_t
ProductDirichlet::remove_observations(const word_t& word) {
        printf("DIST: removing %d observations\n", (int)(word.length/counts->size1));
        for (size_t i = 0; i < word.length; i += counts->size1) {
                for (size_t j = 0; j < counts->size1; j++) {
                        switch ((*word.sequences)[word.sequence][word.position+i+j]) {
                        case 'A':
                        case 'a':
                                gsl_matrix_set(counts, j, 0,
                                               gsl_matrix_get(counts, j, 0)-1);
                        break;
                        case 'C':
                        case 'c':
                                gsl_matrix_set(counts, j, 1,
                                               gsl_matrix_get(counts, j, 1)-1);
                        break;
                        case 'G':
                        case 'g':
                                gsl_matrix_set(counts, j, 2,
                                               gsl_matrix_get(counts, j, 2)-1);
                        break;
                        case 'T':
                        case 't':
                                gsl_matrix_set(counts, j, 3,
                                               gsl_matrix_get(counts, j, 3)-1);
                        break;
                        }
                }
        }
        return word.length/counts->size1;
}

size_t
ProductDirichlet::add_observations(const word_t& word) {
        printf("DIST: adding %d observations\n", (int)(word.length/counts->size1));
        for (size_t i = 0; i < word.length; i += counts->size1) {
                for (size_t j = 0; j < counts->size1; j++) {
                        switch ((*word.sequences)[word.sequence][word.position+i+j]) {
                        case 'A':
                        case 'a':
                                gsl_matrix_set(counts, j, 0,
                                               gsl_matrix_get(counts, j, 0)+1);
                        break;
                        case 'C':
                        case 'c':
                                gsl_matrix_set(counts, j, 1,
                                               gsl_matrix_get(counts, j, 1)+1);
                        break;
                        case 'G':
                        case 'g':
                                gsl_matrix_set(counts, j, 2,
                                               gsl_matrix_get(counts, j, 2)+1);
                        break;
                        case 'T':
                        case 't':
                                gsl_matrix_set(counts, j, 3,
                                               gsl_matrix_get(counts, j, 3)+1);
                        break;
                        }
                }
        }
        return word.length/counts->size1;
}

double ProductDirichlet::pdf(const word_t& word) const {
        double result = 1;

        for (size_t i = 0; i < word.length; i += counts->size1) {
                for (size_t j = 0; j < counts->size1; j++) {
                        double sum = 0;
                        for (size_t k = 0; k < this->counts->size2; k++) {
                                sum += gsl_matrix_get(this->counts, j, k);
                        }
                        switch ((*word.sequences)[word.sequence][word.position+i+j]) {
                        case 'A':
                        case 'a':
                                result *= gsl_matrix_get(this->counts, j, 0)/sum;
                        break;
                        case 'C':
                        case 'c':
                                result *= gsl_matrix_get(this->counts, j, 1)/sum;
                        break;
                        case 'G':
                        case 'g':
                                result *= gsl_matrix_get(this->counts, j, 2)/sum;
                        break;
                        case 'T':
                        case 't':
                                result *= gsl_matrix_get(this->counts, j, 3)/sum;
                        break;
                        }
                }
        }
        return result;
}
