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
#include <distribution.hh>

using namespace std;

gsl_rng* _r;

ProductDirichlet::ProductDirichlet(gsl_matrix* _alpha)
{
        for (size_t i = 0; i < _alpha->size1; i++) {
                double sum = 0;
                alpha.push_back (vector<double>(_alpha->size2+1, 0));
                counts.push_back(vector<size_t>(_alpha->size2+1, 0));
                for (size_t j = 0; j < _alpha->size2; j++) {
                        alpha[i][j] = gsl_matrix_get(_alpha, i, j);
                        sum += alpha[i][j];
                }
                alpha[i][_alpha->size2] = sum;
        }
}

ProductDirichlet::ProductDirichlet(const ProductDirichlet& distribution)
        : alpha(distribution.alpha), counts(distribution.counts)
{
}

ProductDirichlet::~ProductDirichlet() {
}

ProductDirichlet*
ProductDirichlet::clone() const {
        return new ProductDirichlet(*this);
}

size_t
ProductDirichlet::count_observations(const word_t& word) const {
        return word.length/counts.size();
}

size_t
ProductDirichlet::remove_observations(const word_t& word) {
        for (size_t i = 0; i < word.length; i += counts.size()) {
                for (size_t j = 0; j < counts.size(); j++) {
                        const char index = word.sequences[word.sequence][word.position+i+j];
                        counts[j][index]--;
                        counts[j][4]--;
                }
        }
        return word.length/counts.size();
}

size_t
ProductDirichlet::add_observations(const word_t& word) {
        for (size_t i = 0; i < word.length; i += counts.size()) {
                for (size_t j = 0; j < counts.size(); j++) {
                        const char index = word.sequences[word.sequence][word.position+i+j];
                        counts[j][index]++;
                        counts[j][4]++;
                }
        }
        return word.length/counts.size();
}

double ProductDirichlet::pdf(const word_t& word) const {
        double result = 1;

        for (size_t i = 0; i < word.length; i += counts.size()) {
                for (size_t j = 0; j < counts.size(); j++) {
                        const char index = word.sequences[word.sequence][word.position+i+j];
                        result *= (counts[j][index]+alpha[j][index])/(counts[j][4]+alpha[j][4]);
                }
        }
        return result;
}

double ProductDirichlet::log_pdf(const word_t& word) const {
        double result = 0;

        for (size_t i = 0; i < word.length; i += counts.size()) {
                for (size_t j = 0; j < counts.size(); j++) {
                        const char index = word.sequences[word.sequence][word.position+i+j];
                        result += log((counts[j][index]+alpha[j][index])/(counts[j][4]+alpha[j][4]));
                }
        }
        return result;
}

//
// \sum_{x \in X} n_x log(\frac{n_x + \alpha_x}{\sum_{x' \in X} n_{x'} + \alpha_{x'}})
//
double ProductDirichlet::log_likelihood() const {
        double result = 0;

        for (size_t j = 0; j < counts.size(); j++) {
                for (size_t k = 0; k < counts[j].size(); k++) {
                        result += counts[j][k]*log((counts[j][k]+alpha[j][k])/(counts[j][4]+alpha[j][4]));
                }
        }
        return result;
}

////////////////////////////////////////////////////////////////////////////////

BivariateNormal::BivariateNormal(
        const gsl_matrix* Sigma,
        const gsl_matrix* Sigma_0,
        const gsl_vector* mu_0)
        : _N(0), _dimension(2);
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
        _Sigma_inv    = gsl_matrix_alloc(_dimension, _dimension);
        _mu           = gsl_vector_calloc(_dimension);

        // alloc posterior
        _Sigma_N      = gsl_matrix_alloc(_dimension, _dimension);
        _Sigma_N_inv  = gsl_matrix_alloc(_dimension, _dimension);

        // prior
        inverse(_Sigma_0_inv, Sigma_0);
        gsl_vector_memcpy(_mu_0, mu_0);

        // likelihood
        inverse(_Sigma_inv, Sigma);

        // posterior
        gsl_matrix_memcpy(_Sigma_N, Sigma_0);
        inverse(_Sigma_N_inv, _Sigma_N);
}

BivariateNormal::~BivariateNormal()
{
        // tmp
        gsl_matrix_free(_tmp_env);
        gsl_permutation_free(_inv_perm);
        gsl_vector_free(_tmp1);
        gsl_vector_free(_tmp2);

        // prior
        gsl_matrix_free(_Sigma_0_inv);
        gsl_vector_free(_mu_0);

        // likelihood
        gsl_matrix_free(_Sigma_inv);
        gsl_vector_free(_mu);

        // posterior
        gsl_matrix_free(_Sigma_N);
        gsl_matrix_free(_Sigma_N_inv);
        gsl_vector_free(_mu_N);
}

void
BivariateNormal::inverse(gsl_matrix* dst, gsl_matrix* src) {
        int _inv_s = 0;
        gsl_matrix_memcpy(_inv_tmp, src);
        gsl_linalg_LU_decomp(_inv_tmp, _inv_perm, &_inv_s);
        gsl_linalg_LU_invert(_inv_tmp, _inv_perm, dst);
}

void
BivariateNormal::update()
{
        // update _mu_N and _Sigma_N

        // compute posterior covariance _Sigma)N and _Sigma_N_inv
        gsl_matrix_memcpy(_Sigma_N_inv, _Sigma_inv);
        gsl_matrix_scale(_Sigma_N_inv, _N);
        gsl_matrix_add(_Sigma_N_inv, _Sigma_0_inv);
        inverse(_Sigma_N, _Sigma_N_inv);

        // compute posterior mean _mu_N
        // tmp1 = Sigma^-1 mu
        gsl_blas_dgemv(CblasNoTrans, 1.0, _Sigma_inv,   _mu,   0.0, tmp1);
        // tmp2 = Sigma_0^-1 mu_0
        gsl_blas_dgemv(CblasNoTrans, 1.0, _Sigma_inv_0, _mu_0, 0.0, tmp2);
        // tmp1 = N Sigma^-1 mu
        gsl_vector_scale(tmp1, _N);
        // tmp1 = N Sigma^-1 mu + Sigma_0^-1 mu_0
        gsl_vector_add(tmp1, tmp2);
        // _mu_N = Sigma_N (N Sigma^-1 mu + Sigma_0^-1 mu_0)
        gsl_blas_dgemv(CblasTrans, 1.0, _Sigma_N, tmp1, 0.0, _mu_N);

        gsl_vector_free(tmp1);
        gsl_vector_free(tmp2);

}

size_t
BivariateNormal::add_observations(const vector<double>& x)
{
        for (size_t i = 0; i < _dimension; i++) {
                gsl_vector_set(_mu, i, (_N*gsl_vector_get(_mu, i)+x[i])/(N+1))
        }
        _N++;
        update();
}

size_t
BivariateNormal::remove_observations(const vector<double>& x)
{
        for (size_t i = 0; i < _dimension; i++) {
                gsl_vector_set(_mu, i, (_N*gsl_vector_get(_mu, i)-x[i])/(N-1))
        }
        _N--;
        update();
}

size_t
BivariateNormal::count_observations(const vector<double>& x) const
{
        return 1;
}

double BivariateNormal::pdf(const vector<double>& x) const {
        double mu_x = gsl_vector_get(_mu_N, 0);
        double mu_y = gsl_vector_get(_mu_N, 1);
        double sigma_x = sqrt(gsl_matrix_get(_Sigma_N, 0, 0));
        double sigma_y = sqrt(gsl_matrix_get(_Sigma_N, 1, 1));
        double rho  = gsl_matrix_get(_Sigma_N, 0, 1)/(sigma_x*sigma_y);

        return gsl_ran_bivariate_gaussian_pdf(x[0]-mu_x, x[1]-mu_y, sigma_x, sigma_y, rho);
}

double BivariateNormal::log_pdf(const vector<double>& x) const {
        return log(pdf(x));
}

double
BivariateNormal::log_likelihood() const
{
        return 0;
}
