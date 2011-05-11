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

#include <vector>

#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include "data.hh"
#include "statistics.hh"

using namespace std;

gsl_rng* _r;

ProductDirichlet::ProductDirichlet()
{
        update(NULL);
}

ProductDirichlet::ProductDirichlet(
        gsl_matrix* alpha)
{
        update(alpha);
}

ProductDirichlet::~ProductDirichlet() {
        if (this->alpha) {
                gsl_matrix_free(this->alpha);
        }
}

void ProductDirichlet::update(
        gsl_matrix* alpha)
{
        if (this->alpha) {
                gsl_matrix_free(this->alpha);
        }

        this->alpha = alpha;
}

double ProductDirichlet::pdf(char *buf) {
        double result = 1;

        for (unsigned int i = 0; i < this->alpha->size1; i++) {
                double sum = 0;
                for (unsigned int j = 0; j < this->alpha->size2; j++) {
                        sum += gsl_matrix_get(this->alpha, i, j);
                }
                switch (buf[i]) {
                case 'A':
                case 'a':
                        result *= gsl_matrix_get(this->alpha, i, 0)/sum;
                break;
                case 'C':
                case 'c':
                        result *= gsl_matrix_get(this->alpha, i, 1)/sum;
                break;
                case 'G':
                case 'g':
                        result *= gsl_matrix_get(this->alpha, i, 2)/sum;
                break;
                case 'T':
                case 't':
                        result *= gsl_matrix_get(this->alpha, i, 3)/sum;
                break;
                }
        }
        return result;
}
