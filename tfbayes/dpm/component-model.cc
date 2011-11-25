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
#include <string.h>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include <data.hh>
#include <component-model.hh>

using namespace std;

gsl_rng* _r;

ProductDirichlet::ProductDirichlet(gsl_matrix* _alpha, const sequence_data_t<short>& data)
        : _data(data), _size1(_alpha->size1), _size2(_alpha->size2)
{
        for (size_t i = 0; i < _alpha->size1; i++) {
                size_t sum = 0;
                alpha.push_back (vector<size_t>(_alpha->size2+1, 0));
                counts.push_back(vector<size_t>(_alpha->size2+1, 0));
                for (size_t j = 0; j < _alpha->size2; j++) {
                        alpha[i][j] = gsl_matrix_get(_alpha, i, j);
                        sum += alpha[i][j];
                }
                alpha[i][_alpha->size2] = sum;
        }
}

ProductDirichlet::ProductDirichlet(const ProductDirichlet& distribution)
        : alpha(distribution.alpha), counts(distribution.counts), _data(distribution._data),
          _size1(distribution._size1), _size2(distribution._size2)
{
}

ProductDirichlet::~ProductDirichlet() {
}

ProductDirichlet*
ProductDirichlet::clone() const {
        return new ProductDirichlet(*this);
}

size_t
ProductDirichlet::add(const range_t& range) {
        const size_t sequence = range.from[0];
        const size_t from     = range.from[1];
        const size_t to       = range.to[1];
        size_t i;

        for (i = 0; i <= to-from; i++) {
                counts[i%_size1][_data[index_t(sequence, from+i)]]++;
                counts[i%_size1][_size2]++;
        }
        return i/_size1;
}

size_t
ProductDirichlet::remove(const range_t& range) {
        const size_t sequence = range.from[0];
        const size_t from     = range.from[1];
        const size_t to       = range.to[1];
        size_t i;

        for (i = 0; i <= to-from; i++) {
                counts[i%_size1][_data[index_t(sequence, from+i)]]--;
                counts[i%_size1][_size2]--;
        }
        return i/_size1;
}

size_t
ProductDirichlet::count(const range_t& range) {
        return (range.to[1]-range.from[1]+1)/_size1;
}

double ProductDirichlet::predictive(const range_t& range) const {
        const size_t sequence = range.from[0];
        const size_t from     = range.from[1];
        const size_t to       = range.to[1];
        double result = 1;

        for (size_t i = 0; i <= to-from; i++) {
                const char code = _data[index_t(sequence, from+i)];
                result *= (counts[i%_size1][code  ]+alpha[i%_size1][code  ])
                         /(counts[i%_size1][_size2]+alpha[i%_size1][_size2]);
        }
        return result;
}

double ProductDirichlet::log_predictive(const range_t& range) const {
        const size_t sequence = range.from[0];
        const size_t from     = range.from[1];
        const size_t to       = range.to[1];
        double result = 0;

        for (size_t i = 0; i <= to-from; i++) {
                const char code = _data[index_t(sequence, from+i)];
                result += log(counts[i%_size1][code  ]+alpha[i%_size1][  code])
                        - log(counts[i%_size1][_size2]+alpha[i%_size1][_size2]);
        }
        return result;
}

//
// \sum_{x \in X} n_x log(\frac{n_x + \alpha_x}{\sum_{x' \in X} n_{x'} + \alpha_{x'}})
//
double ProductDirichlet::log_likelihood() const {
        double result = 0;

        for (size_t j = 0; j < _size1; j++) {
                for (size_t k = 0; k < _size2; k++) {
                        if (alpha[j][k] != 0) {
                                result += counts[j][k]*(log(counts[j][k]+alpha[j][k]) - log(counts[j][_size2]+alpha[j][_size2]));
                        }
                }
        }
        return result;
}

////////////////////////////////////////////////////////////////////////////////

static
void pt_init(static_pars_tree_t* pt) {
        for (size_t i = 0 ; i < pt->size * pt->as->size ; i++) {
                pt->dirichlet_params[i] = 1.0;
        }
}

ParsimoniousTree::ParsimoniousTree(
        size_t alphabet_size, size_t tree_depth,
        const sequence_data_t<short>& data)
        : _data(data)
{
        _as = as_create(alphabet_size);
        _pt = pt_create(_as, tree_depth);
        _counts_length = pow(alphabet_size, tree_depth+1);
        _counts = (count_t*)malloc(_counts_length*sizeof(count_t));

        /* init data structures */
        pt_init(_pt);
        memset(_counts, 0, _counts_length*sizeof(count_t));
}

ParsimoniousTree::ParsimoniousTree(const ParsimoniousTree& distribution)
        : _counts_length(distribution._counts_length),
          _data(distribution._data)
{
        _as = as_create(distribution._as->size);
        _pt = pt_create(_as, distribution._pt->depth);
        _counts = (count_t*)malloc(_counts_length*sizeof(count_t));

        /* init data structures */
        pt_init(_pt);
        memcpy(_counts, distribution._counts, _counts_length*sizeof(count_t));
}

ParsimoniousTree::~ParsimoniousTree() {
        free(_counts);
        pt_free(_pt);
        as_free(_as);
}

ParsimoniousTree*
ParsimoniousTree::clone() const {
        return new ParsimoniousTree(*this);
}

size_t
ParsimoniousTree::add(const range_t& range) {
        return 1;
}

size_t
ParsimoniousTree::remove(const range_t& range) {
        return 1;
}

size_t
ParsimoniousTree::count(const range_t& range) {
        return 1;
}

double ParsimoniousTree::predictive(const range_t& range) const {
        return 1;
}

double ParsimoniousTree::log_predictive(const range_t& range) const {
        return 1;
}

//
// \sum_{x \in X} n_x log(\frac{n_x + \alpha_x}{\sum_{x' \in X} n_{x'} + \alpha_{x'}})
//
double ParsimoniousTree::log_likelihood() const {
        return 1;
}

////////////////////////////////////////////////////////////////////////////////

BivariateNormal::BivariateNormal(
        const gsl_matrix* Sigma,
        const gsl_matrix* Sigma_0,
        const gsl_vector* mu_0,
        const data_t<vector<double> >& data)
        : _N(0), _dimension(2), _data(data)
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
}

BivariateNormal::BivariateNormal(const BivariateNormal& bn)
        : _N(bn._N), _dimension(bn._dimension), _data(bn._data)
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

BivariateNormal::~BivariateNormal()
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

BivariateNormal*
BivariateNormal::clone() const {
        return new BivariateNormal(*this);
}

void
BivariateNormal::inverse(gsl_matrix* dst, const gsl_matrix* src) {
        int _inv_s = 0;
        gsl_matrix_memcpy(_inv_tmp, src);
        gsl_linalg_LU_decomp(_inv_tmp, _inv_perm, &_inv_s);
        gsl_linalg_LU_invert(_inv_tmp, _inv_perm, dst);
}

void
BivariateNormal::update()
{
        // update _mu_N and _Sigma_N

        // compute posterior covariance _Sigma_N and _Sigma_N_inv
        gsl_matrix_memcpy(_Sigma_N_inv, _Sigma_inv);
        gsl_matrix_scale(_Sigma_N_inv, _N);
        gsl_matrix_add(_Sigma_N_inv, _Sigma_0_inv);
        inverse(_Sigma_N, _Sigma_N_inv);

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
BivariateNormal::add(const range_t& range)
{
        const_iterator_t<vector<double> > iterator = _data[range];

        do {
                for (size_t i = 0; i < _dimension; i++) {
                        gsl_vector_set(_mu, i, (_N*gsl_vector_get(_mu, i)+(*iterator)[i])/(_N+1));
                }
                _N++;
        } while(iterator++);

        update();

        return range.to[0]-range.from[0]+1;
}

size_t
BivariateNormal::remove(const range_t& range)
{
        const_iterator_t<vector<double> > iterator = _data[range];

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

        return range.to[0]-range.from[0]+1;
}

size_t
BivariateNormal::count(const range_t& range) {
        return range.to[0]-range.from[0]+1;
}

double BivariateNormal::predictive(const range_t& range) const {
        const_iterator_t<vector<double> > iterator = _data[range];

        double mu_x = gsl_vector_get(_mu_N, 0);
        double mu_y = gsl_vector_get(_mu_N, 1);
        double sigma_x = sqrt(gsl_matrix_get(_Sigma_N, 0, 0));
        double sigma_y = sqrt(gsl_matrix_get(_Sigma_N, 1, 1));
        double rho  = gsl_matrix_get(_Sigma_N, 0, 1)/(sigma_x*sigma_y);
        double result = 1;

        do {
                result *= gsl_ran_bivariate_gaussian_pdf((*iterator)[0]-mu_x, (*iterator)[1]-mu_y, sigma_x, sigma_y, rho);
        } while (iterator++);

        return result;
}

double BivariateNormal::log_predictive(const range_t& range) const {
        return log(predictive(range));
}

double
BivariateNormal::log_likelihood() const
{
        return 0;
}

const gsl_vector*
BivariateNormal::mean() const
{
        return _mu;
}
