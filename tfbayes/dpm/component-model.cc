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

#include <component-model.hh>

using namespace std;

// Multinomial/Dirichlet Model
////////////////////////////////////////////////////////////////////////////////

ProductDirichlet::ProductDirichlet(const matrix<double>& _alpha, const sequence_data_t<short>& data)
        : _data(data), _size1(_alpha.size()), _size2(_alpha[0].size())
{
        for (size_t i = 0; i < _alpha.size(); i++) {
                double sum = 0;
                alpha.push_back (counts_t(_size2+1, 0.0));
                counts.push_back(counts_t(_size2+1, 0.0));
                for (size_t j = 0; j < _alpha[i].size(); j++) {
                        alpha[i][j] = _alpha[i][j];
                        sum        +=  alpha[i][j];
                }
                alpha[i][_alpha[i].size()] = sum;
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
        const size_t sequence = range.index[0];
        const size_t position = range.index[1];
        const size_t length   = range.length;
        size_t i;

        for (i = 0; i < length; i++) {
                counts[i%_size1][_data[seq_index_t(sequence, position+i)]]++;
                counts[i%_size1][_size2]++;
        }
        return i/_size1;
}

size_t
ProductDirichlet::remove(const range_t& range) {
        const size_t sequence = range.index[0];
        const size_t position = range.index[1];
        const size_t length   = range.length;
        size_t i;

        for (i = 0; i < length; i++) {
                counts[i%_size1][_data[seq_index_t(sequence, position+i)]]--;
                counts[i%_size1][_size2]--;
        }
        return i/_size1;
}

size_t
ProductDirichlet::count(const range_t& range) {
        return range.length/_size1;
}

double ProductDirichlet::predictive(const range_t& range) {
        const size_t sequence = range.index[0];
        const size_t position = range.index[1];
        const size_t length   = range.length;
        double result = 1;

        for (size_t i = 0; i < length; i++) {
                const char code = _data[seq_index_t(sequence, position+i)];
                result *= (counts[i%_size1][ code ]+alpha[i%_size1][ code ])
                         /(counts[i%_size1][_size2]+alpha[i%_size1][_size2]);
        }
        return result;
}

double ProductDirichlet::log_predictive(const range_t& range) {
        const size_t sequence = range.index[0];
        const size_t position = range.index[1];
        const size_t length   = range.length;
        double result = 0;

        for (size_t i = 0; i < length; i++) {
                const char code = _data[seq_index_t(sequence, position+i)];
                result += log(counts[i%_size1][ code ]+alpha[i%_size1][ code ])
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

// Markov Chain Mixture
////////////////////////////////////////////////////////////////////////////////

MarkovChainMixture::MarkovChainMixture(
        size_t alphabet_size, size_t max_context,
        const sequence_data_t<short>& data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        cluster_tag_t cluster_tag)
        : _data(data), _cluster_assignments(cluster_assignments),
          _cluster_tag(cluster_tag), _max_context(max_context),
          _alphabet_size(alphabet_size), _entropy_max(log(_alphabet_size))
{
        _length = context_t::counts_size(alphabet_size, max_context);
        _alpha  = (double*)malloc(_length*sizeof(double));
        _counts = (double*)malloc(_length*sizeof(double));
        _counts_sum = (double*)malloc(_length/_alphabet_size*sizeof(double));
        _entropy    = (double*)malloc(_length/_alphabet_size*sizeof(double));
        _parents    = (int*)malloc(_length*sizeof(int));

        /* for likelihood computations */
        _counts_tmp = (double*)malloc(_length*sizeof(double));

        /* init data structures */
        for (size_t i = 0; i < _length; i++) {
                _alpha[i]  = 0.5;
                _counts[i] = 0.0;
        }
        /* init entropies */
        for (size_t i = 0; i < _length/_alphabet_size; i++) {
                _entropy[i]    = 1;
                _counts_sum[i] = 0; /* also init counts_sum to zero */
        }
        /* init counts_sum */
        for (size_t i = 0; i < _length; i++) {
                _counts_sum[i/_alphabet_size] += _alpha[i] + _counts[i];
        }
        /* for each node compute its parent */
        for (size_t i = 0; i < _length/_alphabet_size; i++) {
                for (size_t j = 0; j < _alphabet_size; j++) {
                        if (i == 0) {
                                _parents[j] = -1;
                        }
                        else {
                                _parents[_alphabet_size*i + j] =
                                        _alphabet_size*((i-1)/_alphabet_size) + j;
                        }
                }
        }
        /* compute context */
        for (size_t i = 0; i < _data.size(); i++) {
                const nucleotide_sequence_t& seq = (const nucleotide_sequence_t&)_data[i];
                _context.push_back(seq_context_t(seq, _max_context, alphabet_size));
        }
}

MarkovChainMixture::MarkovChainMixture(const MarkovChainMixture& distribution)
        : _length(distribution._length),
          _data(distribution._data),
          _context(distribution._context),
          _cluster_assignments(distribution._cluster_assignments),
          _cluster_tag(distribution._cluster_tag),
          _max_context(distribution._max_context),
          _alphabet_size(distribution._alphabet_size),
          _entropy_max(distribution._entropy_max)
{
        _alpha      = (double*)malloc(_length*sizeof(double));
        _counts     = (double*)malloc(_length*sizeof(double));
        _counts_sum = (double*)malloc(_length/_alphabet_size*sizeof(double));
        _entropy    = (double*)malloc(_length/_alphabet_size*sizeof(double));
        _parents    = (int*)malloc(_length*sizeof(int));

        /* for likelihood computations */
        _counts_tmp = (double*)malloc(_length*sizeof(double));

        /* init data structures */
        memcpy(_alpha,   distribution._alpha,   _length*sizeof(double));
        memcpy(_counts,  distribution._counts,  _length*sizeof(double));
        memcpy(_counts_sum, distribution._counts_sum, _length/_alphabet_size*sizeof(double));
        memcpy(_entropy, distribution._counts_sum, _length/_alphabet_size*sizeof(double));
        memcpy(_parents, distribution._parents, _length*sizeof(int));
}

MarkovChainMixture::~MarkovChainMixture() {
        free(_alpha);
        free(_counts);
        free(_parents);
        free(_counts_sum);
        free(_counts_tmp);
}

MarkovChainMixture*
MarkovChainMixture::clone() const {
        return new MarkovChainMixture(*this);
}

size_t
MarkovChainMixture::max_from_context(const range_t& range) const
{
        const size_t sequence = range.index[0];
        const size_t position = range.index[1];
        ssize_t from = position;

        for (size_t i = 0; i < _max_context && from > 0 && _cluster_assignments[seq_index_t(sequence, from-1)] == _cluster_tag; i++) {
                from--;
        }

        return position-from;
}

size_t
MarkovChainMixture::max_to_context(const range_t& range) const
{
        const size_t sequence = range.index[0];
        const size_t position = range.index[1];
        const size_t length   = range.length;
        const ssize_t sequence_length = _data.size(sequence);
        ssize_t to = position+length-1;

        for (size_t i = 0; i < _max_context && to < sequence_length-1 && _cluster_assignments[seq_index_t(sequence, to + 1)] == _cluster_tag; i++) {
                to++;
        }

        return to-position-length+1;
}

void
MarkovChainMixture::update_entropy(int code)
{
        const int from = code - (code%_alphabet_size);
        const int k    = code/_alphabet_size;

        _entropy[k] = 0;
        for (size_t i = from; i < from+_alphabet_size; i++) {
                const double p = (_counts[i]+_alpha[i])/_counts_sum[k];
                _entropy[k] -= p*log(p);
        }
        _entropy[k] /= _entropy_max;
}

size_t
MarkovChainMixture::add(const range_t& range) {
        const size_t sequence     = range.index[0];
        const size_t length       = range.length;
        const size_t from_context = max_from_context(range);
        const size_t   to_context = max_to_context(range);

        for (size_t i = 0; i < length+to_context; i++) {
                const size_t pos   = range.index[1]+i;
                const size_t c_max = min(from_context+i, _max_context);
                const size_t c_min = i < length ? 0 : i-(length-1);
                for (size_t c = c_min; c <= c_max; c++) {
                        const int code = _context[sequence][pos][c];
                        if (code != -1) {
                                _counts[code]++;
                                _counts_sum[code/_alphabet_size]++;
                                update_entropy(code);
                        }
                }
        }

        return range.length;
}

size_t
MarkovChainMixture::remove(const range_t& range) {
        const size_t sequence     = range.index[0];
        const size_t length       = range.length;
        const size_t from_context = max_from_context(range);
        const size_t   to_context = max_to_context(range);

        for (size_t i = 0; i < length+to_context; i++) {
                const size_t pos   = range.index[1]+i;
                const size_t c_max = min(from_context+i, _max_context);
                const size_t c_min = i < length ? 0 : i-(length-1);
                for (size_t c = c_min; c <= c_max; c++) {
                        const int code = _context[sequence][pos][c];
                        if (code != -1) {
                                _counts[code]--;
                                _counts_sum[code/_alphabet_size]--;
                                update_entropy(code);
                        }
                }
        }

        return range.length;
}

size_t
MarkovChainMixture::count(const range_t& range) {
        return range.length;
}

double MarkovChainMixture::predictive(const range_t& range) {
        const size_t sequence     = range.index[0];
        const size_t length       = range.length;
        const size_t from_context = max_from_context(range);
        const size_t   to_context = max_to_context(range);
        double result = 1;

        for (size_t i = 0; i < length+to_context; i++) {
                const size_t pos   = range.index[1]+i;
                const size_t c_max = min(from_context+i, _max_context);
                const size_t c_min = i < length ? 0 : i-(length-1);
                double weight = 0;
                double weight_sum = 0;
                double partial_result = 0;
                for (size_t c = c_min; c <= c_max; c++) {
                        const int code = _context[sequence][pos][c];
                        assert(code != -1);
                        // compute mixture weight
                        if (c == c_max) {
                                weight = 1 - weight_sum;
                        }
                        else {
                                weight = (1 - _entropy[code/_alphabet_size])*(1 - weight_sum);
                        }
                        weight_sum += weight;
                        // compute mixture component
                        partial_result +=
                                weight*(_counts[code]+_alpha[code])/
                                _counts_sum[code/_alphabet_size];
                }
                // save result
                result *= partial_result;
        }

        return result;
}

double MarkovChainMixture::log_predictive(const range_t& range) {
        return log(predictive(range));
}

double MarkovChainMixture::log_likelihood(size_t pos) const
{
        const double n = _counts_tmp[pos];
        double weight = 0;
        double weight_sum = 0;
        double result = 0;
        vector<int> path;

        // compute the path up the tree
        for (int i = pos; i != -1; i = _parents[i]) {
                path.push_back(i);
        }
        // compute likelihood for this node
        for (size_t c = 0; c < path.size(); c++) {
                const int code = path[path.size()-c-1];
                // compute mixture weight
                if (c == path.size()) {
                        weight = 1 - weight_sum;
                }
                else {
                        weight = (1 - _entropy[code/_alphabet_size])*(1 - weight_sum);
                }
                weight_sum += weight;
                // compute mixture component
                result += weight*(_counts[code]+_alpha[code])/
                        _counts_sum[code/_alphabet_size];
        }
        return n*log(result);
}

void MarkovChainMixture::substract_counts(size_t pos) const
{
        const double n = _counts_tmp[pos];

        for (int i = pos; i != -1; i = _parents[i]) {
                _counts_tmp[i] -= n;
        }
}

double MarkovChainMixture::log_likelihood() const {
        memcpy(_counts_tmp, _counts, _length*sizeof(double));
        double result = 0;

        for(size_t i = 0; i <= _max_context; i++) {
                const size_t context = _max_context - i;
                const size_t offset_from = context_t::counts_offset(_alphabet_size, context);
                const size_t offset_to   = context_t::counts_offset(_alphabet_size, context+1);

                // compute likelihood for each node
                for(size_t j = offset_from; j < offset_to; j++) {
                        if (_counts_tmp[j] > 0) {
                                result += log_likelihood(j);
                        }
                }
                // propagate used counts up the tree
                for(size_t j = offset_from; j < offset_to; j++) {
                        if (_counts_tmp[j] > 0) {
                                substract_counts(j);
                        }
                }
        }

        return result;
}

// Variable Order Markov Chain
////////////////////////////////////////////////////////////////////////////////

ParsimoniousTree::ParsimoniousTree(
        size_t alphabet_size, size_t tree_depth,
        const sequence_data_t<short>& data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        cluster_tag_t cluster_tag)
        : _data(data), _cluster_assignments(cluster_assignments),
          _cluster_tag(cluster_tag), _tree_depth(tree_depth)
{
        _as = as_create(alphabet_size);
        _pt = pt_create(_as, tree_depth);
        _counts_length = context_t::counts_size(alphabet_size, tree_depth);
        _counts = (count_t*)malloc(_counts_length*sizeof(count_t));

        /* init data structures */
        pt_init(_pt);
        memset(_counts, 0, _counts_length*sizeof(count_t));

        /* compute context */
        for (size_t i = 0; i < _data.size(); i++) {
                const nucleotide_sequence_t& seq = (const nucleotide_sequence_t&)_data[i];
                _context.push_back(seq_context_t(seq, tree_depth, alphabet_size));
        }
}

ParsimoniousTree::ParsimoniousTree(const ParsimoniousTree& distribution)
        : _counts_length(distribution._counts_length),
          _data(distribution._data),
          _context(distribution._context),
          _cluster_assignments(distribution._cluster_assignments),
          _cluster_tag(distribution._cluster_tag),
          _tree_depth(distribution._tree_depth)
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
ParsimoniousTree::max_from_context(const range_t& range) const
{
        const size_t sequence = range.index[0];
        const size_t position = range.index[1];
        ssize_t from = position;

        for (size_t i = 0; i < _tree_depth && from > 0 && _cluster_assignments[seq_index_t(sequence, from-1)] == _cluster_tag; i++) {
                from--;
        }

        return position-from;
}

size_t
ParsimoniousTree::max_to_context(const range_t& range) const
{
        const size_t sequence = range.index[0];
        const size_t position = range.index[1];
        const size_t length   = range.length;
        const ssize_t sequence_length = _data.size(sequence);
        ssize_t to = position+length-1;

        for (size_t i = 0; i < _tree_depth && to < sequence_length-1 && _cluster_assignments[seq_index_t(sequence, to + 1)] == _cluster_tag; i++) {
                to++;
        }

        return to-position-length+1;
}

size_t
ParsimoniousTree::add(const range_t& range) {
        const size_t sequence     = range.index[0];
        const size_t length       = range.length;
        const size_t from_context = max_from_context(range);
        const size_t   to_context = max_to_context(range);

        for (size_t i = 0; i < length+to_context; i++) {
                const size_t pos   = range.index[1]+i;
                const size_t c_max = min(from_context+i, _tree_depth);
                const size_t c_min = i < length ? 0 : i-(length-1);
                for (size_t c = c_min; c <= c_max; c++) {
                        if (_context[sequence][pos][c] != -1) {
                                _counts[_context[sequence][pos][c]]++;
                        }
                }
        }

        return range.length;
}

size_t
ParsimoniousTree::remove(const range_t& range) {
        const size_t sequence     = range.index[0];
        const size_t length       = range.length;
        const size_t from_context = max_from_context(range);
        const size_t   to_context = max_to_context(range);

        for (size_t i = 0; i < length+to_context; i++) {
                const size_t pos   = range.index[1]+i;
                const size_t c_max = min(from_context+i, _tree_depth);
                const size_t c_min = i < length ? 0 : i-(length-1);
                for (size_t c = c_min; c <= c_max; c++) {
                        if (_context[sequence][pos][c] != -1) {
                                _counts[_context[sequence][pos][c]]--;
                        }
                }
        }

        return range.length;
}

size_t
ParsimoniousTree::count(const range_t& range) {
        return range.length;
}

double ParsimoniousTree::predictive(const range_t& range) {
        const double ml1 = pt_ln_marginal_likelihood(_pt, _counts); add(range);
        const double ml2 = pt_ln_marginal_likelihood(_pt, _counts); remove(range);
        return exp(ml2-ml1);
}

double ParsimoniousTree::log_predictive(const range_t& range) {
        const double ml1 = pt_ln_marginal_likelihood(_pt, _counts); add(range);
        const double ml2 = pt_ln_marginal_likelihood(_pt, _counts); remove(range);
        return ml2-ml1;
}

double ParsimoniousTree::log_likelihood() const {
        return 0;
}

// Bivariate Gaussian
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

        return range.length;
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

        return range.length;
}

size_t
BivariateNormal::count(const range_t& range) {
        return range.length;
}

double BivariateNormal::predictive(const range_t& range) {
        const_iterator_t<vector<double> > iterator = _data[range];

        double mu_x = gsl_vector_get(_mu_N, 0);
        double mu_y = gsl_vector_get(_mu_N, 1);
        double sigma_x = sqrt(gsl_matrix_get(_Sigma_N, 0, 0));
        double sigma_y = sqrt(gsl_matrix_get(_Sigma_N, 1, 1));
        double rho = gsl_matrix_get(_Sigma_N, 0, 1)/(sigma_x*sigma_y);
        double result = 1;

        do {
                result *= gsl_ran_bivariate_gaussian_pdf((*iterator)[0]-mu_x, (*iterator)[1]-mu_y, sigma_x, sigma_y, rho);
        } while (iterator++);

        return result;
}

double BivariateNormal::log_predictive(const range_t& range) {
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
