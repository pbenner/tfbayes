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

#include <math.h>

#include <iterator>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <data-gaussian.hh>

using namespace std;

DataGaussian::DataGaussian(size_t cluster, size_t samples, gsl_matrix* Sigma, const double* pi)
        : data_t<std::vector<double> >(generate_samples(cluster, samples, Sigma, pi)),
          _elements(samples), _length(1), _cluster(cluster)
{
}

DataGaussian::~DataGaussian()
{
        for (size_t i = 0; i < _cluster; i++) {
                gsl_vector_free(_mu[i]);
        }
        free(_mu);
}

vector<vector<double> >
DataGaussian::generate_samples(
        size_t cluster,
        size_t samples,
        gsl_matrix* Sigma,
        const double* pi)
{
        gsl_ran_discrete_t* gdd = gsl_ran_discrete_preproc(cluster, pi);
        double sigma_x = sqrt(gsl_matrix_get(Sigma, 0, 0));
        double sigma_y = sqrt(gsl_matrix_get(Sigma, 1, 1));
        double rho  = gsl_matrix_get(Sigma, 0, 1)/(sigma_x*sigma_y);

        // alloc memory for means
        _mu = (gsl_vector**)malloc(cluster*sizeof(gsl_matrix*));
        // final samples
        vector<vector<double> > data;

        // generate a random mean for each cluster
        for (size_t i = 0; i < cluster; i++) {
                _mu[i] = gsl_vector_alloc(2);
                gsl_vector_set(_mu[i], 0, (double)rand()/RAND_MAX);
                gsl_vector_set(_mu[i], 1, (double)rand()/RAND_MAX);
        }

        // generate samples
        double x0;
        double x1;
        for (size_t i = 0; i < samples; i++) {
                vector<double> x(2,0);
                size_t j = gsl_ran_discrete(_r, gdd);
                gsl_ran_bivariate_gaussian(_r, sigma_x, sigma_y, rho, &x0, &x1);
                x[0] = gsl_vector_get(_mu[j], 0) + x0;
                x[1] = gsl_vector_get(_mu[j], 1) + x1;
                data.push_back(x);
                indices.push_back(index_t(i));
        }

        // generate a randomized list of indices
        for (DataGaussian::iterator it = begin(); it != end(); it++) {
                indices_randomized.push_back(&(*it));
        }
        shuffle();

        gsl_ran_discrete_free(gdd);

        return data;
}

void
DataGaussian::shuffle() {
        random_shuffle(indices_randomized.begin(), indices_randomized.end());
}

const index_t&
DataGaussian::operator[](size_t i) const {
        return indices[i];
}

size_t
DataGaussian::elements() const {
        return _elements;
}

size_t
DataGaussian::length() const {
        return _length;
}
