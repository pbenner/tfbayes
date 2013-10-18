/* Copyright (C) 2011-2013 Philipp Benner
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

#include <iterator>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <tfbayes/dpm/data-gaussian.hh>

using namespace std;

data_gaussian_t::data_gaussian_t(
        size_t samples,
        const matrix<double>& Sigma,
        const vector<double>& pi)
        : data_t<std::vector<double> >(),
          _elements(samples),
          _length  (1),
          _cluster (pi.size())
{
        gsl_rng* _r = gsl_rng_alloc (gsl_rng_default);

        /* copy contents of pi into a standard array */
        double _pi[pi.size()];
        copy(pi.begin(), pi.end(), _pi);

        /* initialize gsl random generator */
        gsl_ran_discrete_t* gdd = gsl_ran_discrete_preproc(pi.size(), _pi);

        /* covariance matrix for normal samples */
        double sigma_x = sqrt(Sigma[0][0]);
        double sigma_y = sqrt(Sigma[1][1]);
        double rho     = Sigma[0][1]/(sigma_x*sigma_y);

        /* final samples */
        matrix<double> data;

        /* temporary storage */
        size_t j;
        double x0;
        double x1;
        vector<double> x(2,0);

        /* generate a random mean for each cluster */
        _mu = matrix<double>(pi.size(), 2);
        for (size_t i = 0; i < pi.size(); i++) {
                _mu[i][0] = 1.2*(double)rand()/RAND_MAX-0.5;
                _mu[i][1] = 1.2*(double)rand()/RAND_MAX-0.5;
        }
        /* generate samples */
        for (size_t i = 0; i < samples; i++) {
                /* select component */
                j = gsl_ran_discrete(_r, gdd);
                /* normal sample */
                gsl_ran_bivariate_gaussian(_r, sigma_x, sigma_y, rho, &x0, &x1);
                x[0] = _mu[j][0] + x0;
                x[1] = _mu[j][1] + x1;
                /* save data point */
                data.push_back(x);
                /* save cluster assignment */
                _initial_cluster_assignments.push_back(j);
        }
        gsl_ran_discrete_free(gdd);

        /* copy data */
        data_t<std::vector<double> >::operator=(data);

        /* initialize indices */
        for (size_t i = 0; i < samples; i++) {
                index_t* index = new index_t(i);
                indices.push_back(index);
                sampling_indices.push_back(index);
        }
        shuffle();

        /* free random number generator */
        gsl_rng_free (_r);
}

data_gaussian_t::data_gaussian_t(const data_gaussian_t& data)
        : _elements(data._elements),
          _length  (data._length),
          _cluster (data._cluster),
          _mu      (data._mu),
          _initial_cluster_assignments(data._initial_cluster_assignments)
{
        // generate a randomized list of indices
        for (data_gaussian_t::const_iterator it = data.begin(); it != data.end(); it++) {
                index_i* index = (**it).clone();
                indices.push_back(index);
                sampling_indices.push_back(index);
        }
        shuffle();
}

data_gaussian_t::~data_gaussian_t()
{
        for (iterator it = begin(); it != end(); it++) {
                delete(*it);
        }
}

void
data_gaussian_t::shuffle() {
        random_shuffle(sampling_indices.begin(), sampling_indices.end());
}

size_t
data_gaussian_t::elements() const {
        return _elements;
}

size_t
data_gaussian_t::length() const {
        return _length;
}

const matrix<double>&
data_gaussian_t::initial_means() const {
        return _mu;
}

const vector<double>&
data_gaussian_t::initial_cluster_assignments() const {
        return _initial_cluster_assignments;
}
