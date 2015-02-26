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

#include <cmath>        // std::exp, std::log
#include <string>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/fastarithmetics/fast-lnbeta.hh>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

bivariate_normal_t::bivariate_normal_t(
        const std::matrix<double>& arg_Sigma,
        const std::matrix<double>& arg_Sigma_0,
        const std::vector<double>& arg_mu_0,
        const data_t<vector<double> >& data)
        : component_model_t(),
          _N               (0),
          _dimension       (2),
          _data            (&data)
{
        gsl_matrix* Sigma   = to_gsl_matrix(arg_Sigma);
        gsl_matrix* Sigma_0 = to_gsl_matrix(arg_Sigma_0);
        gsl_vector* mu_0    = to_gsl_vector(arg_mu_0);

        // alloc tmp
        _inv_tmp  = gsl_matrix_alloc(_dimension, _dimension);
        _inv_perm = gsl_permutation_alloc(_dimension);
        _tmp1     = gsl_vector_alloc(_dimension);
        _tmp2     = gsl_vector_alloc(_dimension);
 
        // alloc prior
        _Sigma_0_inv  = gsl_matrix_alloc(_dimension, _dimension);
        _mu_0         = gsl_vector_alloc(_dimension);

        // alloc likelihood
        _Sigma        = gsl_matrix_alloc(_dimension, _dimension);
        _Sigma_inv    = gsl_matrix_alloc(_dimension, _dimension);
        _mu           = gsl_vector_calloc(_dimension);

        // alloc posterior
        _Sigma_N      = gsl_matrix_alloc(_dimension, _dimension);
        _Sigma_N_inv  = gsl_matrix_alloc(_dimension, _dimension);
        _mu_N         = gsl_vector_alloc(_dimension);
 
        // prior
        inverse(_Sigma_0_inv, Sigma_0);
        gsl_vector_memcpy(_mu_0, mu_0);

        // likelihood
        gsl_matrix_memcpy(_Sigma, Sigma);
        inverse(_Sigma_inv, Sigma);

        // posterior
        gsl_matrix_memcpy(_Sigma_N, Sigma_0);
        gsl_matrix_add   (_Sigma_N, Sigma);
        inverse(_Sigma_N_inv, _Sigma_N);
        gsl_vector_memcpy(_mu_N, mu_0);

        // free args
        gsl_matrix_free(Sigma);
        gsl_matrix_free(Sigma_0);
        gsl_vector_free(mu_0);
}

bivariate_normal_t::bivariate_normal_t(const bivariate_normal_t& bn)
        : component_model_t(bn),
          _N               (bn._N),
          _dimension       (bn._dimension),
          _data            (bn._data)
{
        // alloc tmp
        _inv_tmp  = gsl_matrix_alloc(_dimension, _dimension);
        _inv_perm = gsl_permutation_alloc(_dimension);
        _tmp1     = gsl_vector_alloc(_dimension);
        _tmp2     = gsl_vector_alloc(_dimension);

        // alloc prior
        _Sigma_0_inv  = gsl_matrix_alloc(_dimension, _dimension);
        _mu_0         = gsl_vector_alloc(_dimension);

        // alloc likelihood
        _Sigma        = gsl_matrix_alloc(_dimension, _dimension);
        _Sigma_inv    = gsl_matrix_alloc(_dimension, _dimension);
        _mu           = gsl_vector_calloc(_dimension);

        // alloc posterior
        _mu_N         = gsl_vector_alloc(_dimension);
        _Sigma_N      = gsl_matrix_alloc(_dimension, _dimension);
        _Sigma_N_inv  = gsl_matrix_alloc(_dimension, _dimension);

        // prior
        gsl_matrix_memcpy(_Sigma_0_inv, bn._Sigma_0_inv);
        gsl_vector_memcpy(_mu_0,        bn._mu_0);

        // likelihood
        gsl_matrix_memcpy(_Sigma,       bn._Sigma);
        gsl_matrix_memcpy(_Sigma_inv,   bn._Sigma_inv);
        gsl_vector_memcpy(_mu,          bn._mu);

        // posterior
        gsl_matrix_memcpy(_Sigma_N,     bn._Sigma_N);
        gsl_matrix_memcpy(_Sigma_N_inv, bn._Sigma_N_inv);
        gsl_vector_memcpy(_mu_N,        bn._mu_N);
}

bivariate_normal_t::~bivariate_normal_t()
{
        // tmp
        gsl_matrix_free(_inv_tmp);
        gsl_permutation_free(_inv_perm);
        gsl_vector_free(_tmp1);
        gsl_vector_free(_tmp2);

        // prior
        gsl_matrix_free(_Sigma_0_inv);
        gsl_vector_free(_mu_0);

        // likelihood
        gsl_matrix_free(_Sigma);
        gsl_matrix_free(_Sigma_inv);
        gsl_vector_free(_mu);

        // posterior
        gsl_matrix_free(_Sigma_N);
        gsl_matrix_free(_Sigma_N_inv);
        gsl_vector_free(_mu_N);
}

bivariate_normal_t*
bivariate_normal_t::clone() const {
        return new bivariate_normal_t(*this);
}

bivariate_normal_t&
bivariate_normal_t::operator=(const component_model_t& component_model)
{
        bivariate_normal_t tmp(static_cast<const bivariate_normal_t&>(component_model));
        swap(*this, tmp);
        return *this;
}

void
bivariate_normal_t::inverse(gsl_matrix* dst, const gsl_matrix* src) {
        int _inv_s = 0;
        gsl_matrix_memcpy(_inv_tmp, src);
        gsl_linalg_LU_decomp(_inv_tmp, _inv_perm, &_inv_s);
        gsl_linalg_LU_invert(_inv_tmp, _inv_perm, dst);
}

void
bivariate_normal_t::update()
{
        // update _mu_N and _Sigma_N

        // compute posterior covariance _Sigma_N and _Sigma_N_inv
        gsl_matrix_memcpy(_Sigma_N_inv, _Sigma_inv);
        gsl_matrix_scale (_Sigma_N_inv, _N);
        gsl_matrix_add   (_Sigma_N_inv, _Sigma_0_inv);
        inverse          (_Sigma_N, _Sigma_N_inv);

        // compute posterior mean _mu_N
        // tmp1 = Sigma^-1 mu
        gsl_blas_dgemv(CblasNoTrans, 1.0, _Sigma_inv,   _mu,   0.0, _tmp1);
        // tmp2 = Sigma_0^-1 mu_0
        gsl_blas_dgemv(CblasNoTrans, 1.0, _Sigma_0_inv, _mu_0, 0.0, _tmp2);
        // tmp1 = N Sigma^-1 mu
        gsl_vector_scale(_tmp1, _N);
        // tmp1 = N Sigma^-1 mu + Sigma_0^-1 mu_0
        gsl_vector_add(_tmp1, _tmp2);
        // _mu_N = Sigma_N (N Sigma^-1 mu + Sigma_0^-1 mu_0)
        gsl_blas_dgemv(CblasTrans, 1.0, _Sigma_N, _tmp1, 0.0, _mu_N);

        gsl_matrix_add(_Sigma_N, _Sigma);
}

size_t
bivariate_normal_t::add(const range_t& range)
{
        const_iterator_t<vector<double> > iterator = data()[range];

        do {
                for (size_t i = 0; i < _dimension; i++) {
                        gsl_vector_set(_mu, i, (_N*gsl_vector_get(_mu, i)+(*iterator)[i])/(_N+1));
                }
                _N++;
        } while(iterator++);

        update();

        return range.length();
}

size_t
bivariate_normal_t::remove(const range_t& range)
{
        const_iterator_t<vector<double> > iterator = data()[range];

        do {
                for (size_t i = 0; i < _dimension; i++) {
                        if (_N - 1 == 0) {
                                gsl_vector_set(_mu, i, 0);
                        }
                        else {
                                gsl_vector_set(_mu, i, (_N*gsl_vector_get(_mu, i)-(*iterator)[i])/(_N-1));
                        }
                }
                _N--;
        } while(iterator++);

        update();

        return range.length();
}

size_t
bivariate_normal_t::count(const range_t& range) {
        return range.length();
}

double bivariate_normal_t::predictive(const range_t& range) {
        const_iterator_t<vector<double> > iterator = data()[range];

        double mu_x    = gsl_vector_get(_mu_N, 0);
        double mu_y    = gsl_vector_get(_mu_N, 1);
        double sigma_x = sqrt(gsl_matrix_get(_Sigma_N, 0, 0));
        double sigma_y = sqrt(gsl_matrix_get(_Sigma_N, 1, 1));
        double rho     = gsl_matrix_get(_Sigma_N, 0, 1)/(sigma_x*sigma_y);
        double result  = 1;

        do {
                result *= gsl_ran_bivariate_gaussian_pdf((*iterator)[0]-mu_x, (*iterator)[1]-mu_y, sigma_x, sigma_y, rho);
        } while (iterator++);

        return result;
}

double bivariate_normal_t::predictive(const vector<range_t>& range_set) {
        assert(false);
}

double bivariate_normal_t::log_predictive(const range_t& range) {
        return log(predictive(range));
}

double bivariate_normal_t::log_predictive(const vector<range_t>& range_set) {
        return log(predictive(range_set));
}

double
bivariate_normal_t::log_likelihood() const
{
        return 0;
}

std::vector<double>
bivariate_normal_t::mean() const
{
        return from_gsl_vector(_mu);
}
