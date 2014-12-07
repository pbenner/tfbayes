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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>    // std::count
#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_psi.h>

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/fastarithmetics/fast-lnbeta.hh>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/unordered_map.hpp> 

using namespace std;

// Methods to integrate out the pseudocounts of a Dirichlet
// distribution, where we use a Gamma hyperprior
////////////////////////////////////////////////////////////////////////////////

struct gamma_marginal_data {
        const independence_background_t::counts_t& counts;
        const independence_background_t::counts_t& alpha;
        const boost::math::gamma_distribution<> distribution;
};

double
gamma_marginal_f(double * x, size_t dim, void * params)
{
        /* store the result on normal scale */
        double result;
        /* casted parameters */
        struct gamma_marginal_data* data = (gamma_marginal_data *)params;

        /* lnbeta(n, a) - lnbeta(a) */
        result = exp(fast_lnbeta<data_tfbs_t::alphabet_size>(data->counts, x) -
                     fast_lnbeta<data_tfbs_t::alphabet_size>(x));

        /* multiply with gamma distribution */
        for (size_t i = 0; i < data_tfbs_t::alphabet_size; i++) {
                if (data->alpha[i] == -1) {
                        result *= boost::math::pdf(data->distribution, x[i]);
                }
        }
        return result;
}

/* use a gamma prior distribution for the Dirichlet pseudocounts to
 * integrate them out (there is no closed form solution so we need to
 * do this numerically) */
double
gamma_marginal(
        const independence_background_t::counts_t& counts,
        const independence_background_t::counts_t& alpha,
        const double k, const double g)
{
        size_t dim = count(alpha.begin(), alpha.end(), -1);
        double xl[data_tfbs_t::alphabet_size];
        double xu[data_tfbs_t::alphabet_size];
        const gsl_rng_type *T;
        gsl_rng *r;

        size_t calls = 500000;
        double result, err;

        gsl_monte_function F;

        struct gamma_marginal_data data = {
                counts, alpha, boost::math::gamma_distribution<>(k, g)
        };

        for (size_t i = 0; i < dim; i++) {
                        xl[i] =   0.0;
                        xu[i] = 100.0;
        }

        F.f      = gamma_marginal_f;
        F.dim    = dim;
        F.params = &data;
     
        gsl_rng_env_setup();
     
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);

        gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim);
        gsl_monte_miser_integrate(&F, xl, xu, dim, calls, r, s,
                                  &result, &err);
        gsl_monte_miser_free(s);

        return log(result);
}

namespace boost {
        // add hash_value for boost arrays
        static inline
        size_t hash_value(const data_tfbs_t::code_t& counts)
        {
                boost::hash<double> hasher;
                double result = 0.0;

                for (size_t i = 0; i < data_tfbs_t::alphabet_size; i++) {
                        result += hasher(counts[i]);
                }

                return result;
        }

        // overwrite equality operator
        bool operator==(const data_tfbs_t::code_t& counts1, const data_tfbs_t::code_t& _counts2)
        {
                data_tfbs_t::code_t counts2(_counts2);

                for (size_t i = 0; i < data_tfbs_t::alphabet_size; i++) {
                        bool permutation = false;
                        for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                                if (counts1[i] == counts2[j]) {
                                        // yeah, this value exists
                                        permutation = true;
                                        // make sure this value isn't used twice
                                        counts2[i]  = -1.0;
                                }
                        }
                        if (!permutation) {
                                return false;
                        }
                }
                return true;
        }
}

double
hashed_gamma_marginal(
        const independence_background_t::counts_t& counts,
        const independence_background_t::counts_t& alpha,
        const double k, const double g,
        boost::unordered_map<data_tfbs_t::code_t, double>& map,
        boost::shared_mutex& mutex)
{
        {
                // get read access to the map
                boost::shared_lock<boost::shared_mutex> lock(mutex);
                boost::unordered_map<data_tfbs_t::code_t, double>::iterator it = map.find(counts);

                if (it != map.end()) {
                        return it->second;
                }
                // release lock
        }
        {
                // compute value
                double result = gamma_marginal(counts, alpha, k, g);
                // get write access
                boost::unique_lock<boost::shared_mutex> lock(mutex);
                map[counts] = result;

                return result;
        }
}

class background_cache_t {
private:
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
                ar & alpha;
                ar & k;
                ar & g;
                ar & data;
                ar & precomputed_marginal;
        }
public:
        background_cache_t()
                { }
        background_cache_t(
                const independence_background_t::counts_t& alpha,
                const double k, const double g,
                const sequence_data_t<data_tfbs_t::code_t>& data,
                const sequence_data_t<double>& precomputed_marginal)
                : alpha(alpha), k(k), g(g), data(data),
                  precomputed_marginal(precomputed_marginal)
                { }
        bool consistent(
                const independence_background_t::counts_t& alpha,
                const double k, const double g,
                const sequence_data_t<data_tfbs_t::code_t>& data) {
                if (this->alpha.size() != alpha.size()) {
                        return false;
                }
                for (size_t i = 0; i < alpha.size(); i++) {
                        if (this->alpha[i] != alpha[i]) {
                                return false;
                        }
                }
                if (this->k != k || this->g != g) {
                        return false;
                }
                if (this->data.size() != data.size()) {
                        return false;
                }
                for (size_t i = 0; i < data.size(); i++) {
                        if (this->data[i].size() != data[i].size()) {
                                return false;
                        }
                        for (size_t j = 0; j < data[i].size(); j++) {
                                if (this->data[i][j] != data[i][j]) {
                                        return false;
                                }
                        }
                }
                return true;
        }

        independence_background_t::counts_t alpha;
        double k;
        double g;
        sequence_data_t<data_tfbs_t::code_t> data;
        sequence_data_t<double> precomputed_marginal;
};

// Independence Background Model
////////////////////////////////////////////////////////////////////////////////

independence_background_t::independence_background_t(
        const matrix<double>& _alpha,
        const sequence_data_t<data_tfbs_t::code_t>& _data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        boost::optional<const alignment_set_t<>&> alignment_set)
        : component_model_t(cluster_assignments),
          _size(data_tfbs_t::alphabet_size),
          _bg_cluster_tag(0),
          _precomputed_marginal(_data.sizes(), 0),
          _data(&_data)
{
        data_tfbs_t::code_t alpha;

        assert(_alpha.size() == 1);
        assert(_alpha[0].size() == data_tfbs_t::alphabet_size);

        for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                alpha[j] = _alpha[0][j];
        }

        precompute_marginal(alpha);
}

/* This is an independence background model with Dirichlet
 * prior and Gamma distributed pseudocounts. The pseudocounts
 * are numerically integrated out. */
independence_background_t::independence_background_t(
        const matrix<double>& _alpha,
        const double k, const double g,
        const sequence_data_t<data_tfbs_t::code_t>& _data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        thread_pool_t& thread_pool,
        const string& cachefile,
        boost::optional<const alignment_set_t<>&> alignment_set)
        : component_model_t(cluster_assignments),
          _size(data_tfbs_t::alphabet_size),
          _bg_cluster_tag(0),
          _precomputed_marginal(_data.sizes(), 0),
          _data(&_data)
{
        counts_t alpha;

        assert(_alpha.size() == 1);
        assert(_alpha[0].size() == data_tfbs_t::alphabet_size);

        for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                alpha[j] = _alpha[0][j];
        }

        if (std::find(alpha.begin(), alpha.end(), -1) == alpha.end()) {
                // all pseudocounts are given
                precompute_marginal(alpha);
        }
        else {
                // at least one pseudocount has to be integrated
                if (!load_marginal_gamma(alpha, k, g, cachefile)) {
                        precompute_marginal_gamma(alpha, k, g, thread_pool);
                        save_marginal_gamma(alpha, k, g, cachefile);
                }
        }
}

independence_background_t::independence_background_t(const independence_background_t& distribution)
        : component_model_t(distribution),
          _size(distribution._size),
          _bg_cluster_tag(distribution._bg_cluster_tag),
          _precomputed_marginal(distribution._precomputed_marginal),
          _data(distribution._data)
{
}

independence_background_t::~independence_background_t() {
}

independence_background_t*
independence_background_t::clone() const {
        return new independence_background_t(*this);
}

independence_background_t&
independence_background_t::operator=(const component_model_t& component_model)
{
        independence_background_t tmp(
                static_cast<const independence_background_t&>(component_model));
        swap(*this, tmp);
        return *this;
}

bool
independence_background_t::load_marginal_gamma(
        const counts_t& alpha,
        const double k, const double g,
        const string& cachefile)
{
        std::ifstream ifs(cachefile, std::ios::binary);
        background_cache_t background_cache;

        // try to receive precomputed marginals from cache
        if (ifs) {
                boost::archive::binary_iarchive ia(ifs);
                ia >> background_cache;
                if (background_cache.consistent(alpha, k, g, data())) {
                        _precomputed_marginal = background_cache.precomputed_marginal;
                        return true;
                }
                cerr << "Background cache is inconsistent... recomputing!"
                     << endl;
        }
        return false;
}

bool
independence_background_t::save_marginal_gamma(
        const counts_t& alpha,
        const double k, const double g,
        const string& cachefile)
{
        std::ofstream ofs(cachefile);
        if (ofs) {
                boost::archive::binary_oarchive oa(ofs);
                const background_cache_t tmp(alpha, k, g, data(), _precomputed_marginal);
                oa << tmp;
                flockfile(stderr);
                cerr << boost::format("Background cache saved to `%s'.") % cachefile
                     << endl;
                funlockfile(stderr);

                return true;
        }
        return false;
}

void
independence_background_t::precompute_marginal(
        const counts_t& alpha)
{
        /* go through the data and precompute
         * lnbeta(n + alpha) - lnbeta(alpha) */
        for(size_t i = 0; i < data().size(); i++) {
                for(size_t j = 0; j < data()[i].size(); j++) {
                        _precomputed_marginal[i][j] =
                                  fast_lnbeta<data_tfbs_t::alphabet_size>(alpha, data()[i][j])
                                - fast_lnbeta<data_tfbs_t::alphabet_size>(alpha);
                }
        }
}

void
independence_background_t::precompute_marginal_gamma(
        const counts_t& alpha,
        const double k, const double g,
        thread_pool_t& thread_pool)
{
        typedef double (*hgm)(const counts_t&, const counts_t&,
                              const double, const double,
                              boost::unordered_map<counts_t, double>&,
                              boost::shared_mutex&);
        boost::unordered_map<counts_t, double> map;
        boost::shared_mutex mutex;

        flockfile(stderr);
        cerr << "Background gamma shape: " << k << endl
             << "Background gamma scale: " << g << endl
             << endl;
        funlockfile(stderr);
        /* go through the data and precompute
         * lnbeta(n + alpha) - lnbeta(alpha) */
        for(size_t i = 0; i < data().size(); i++) {
                future_vector_t<double> futures(data()[i].size());

                for(size_t j = 0; j < data()[i].size(); j++) {
                        boost::function<double ()> f = boost::bind(
                                static_cast<hgm>(&hashed_gamma_marginal),
                                boost::cref(data()[i][j]), boost::cref(alpha), k, g,
                                boost::ref(map), boost::ref(mutex));

                        futures[j] = thread_pool.schedule(f);
                }
                for(size_t j = 0; j < data()[i].size(); j++) {
                        /* compute percentage by linearly
                         * interpolating two values */
                        const double p = j/(double)data()[i].size();
                        const double q = (p*(i+1.0) + (1.0-p)*i)/(double)data().size();

                        flockfile(stderr);
                        cerr.precision(2);
                        cerr << "\rPrecomputing background... " << setw(6) << fixed
                             << q*100.0 << "%"                  << flush;
                        funlockfile(stderr);

                        _precomputed_marginal[i][j] = futures[j].get();
                }
        }
        flockfile(stderr);
        cout << "\rPrecomputing background...   done." << endl << flush;
        funlockfile(stderr);
}

size_t
independence_background_t::add(const range_t& range) {
        return range.length();
}

size_t
independence_background_t::remove(const range_t& range) {
        return range.length();
}

size_t
independence_background_t::count(const range_t& range) {
        return range.length();
}

/*
 *  p(y|x) = Beta(n(x) + n(y) + alpha) / Beta(n(x) + alpha)
 */
double independence_background_t::predictive(const range_t& range) {
        return exp(log_predictive(range));
}

double independence_background_t::predictive(const vector<range_t>& range_set) {
        return exp(log_predictive(range_set));
}

double independence_background_t::log_predictive(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();
        double result = 0;

        for (size_t i = 0; i < length; i++) {
                const seq_index_t index(sequence, position+i);

                /* counts contains the data count statistic
                 * and the pseudo counts alpha */
                result += _precomputed_marginal[index];
        }

        return result;
}

double independence_background_t::log_predictive(const vector<range_t>& range_set) {
        assert(range_set.size() > 0);

        const size_t length = range_set[0].length();
        double result = 0;

        for (size_t i = 0; i < length; i++) {
                /* loop through all ranges */
                for (size_t k = 0; k < range_set.size(); k++) {

                        const size_t sequence = range_set[k].index()[0];
                        const size_t position = range_set[k].index()[1];
                        const seq_index_t index(sequence, position+i);

                        /* all positions in the alignment are fully
                         * independent, hence we do not need to sum
                         * any counts */
                        result += _precomputed_marginal[index];
                }
        }

        return result;
}

/*
 *  p(x) = Beta(n(x) + alpha) / Beta(alpha)
 */
double independence_background_t::log_likelihood() const {
        double result = 0;

        /* counts contains the data count statistic
         * and the pseudo counts alpha */
        for(size_t i = 0; i < cluster_assignments().size(); i++) {
                for(size_t j = 0; j < cluster_assignments()[i].size(); j++) {
                        if (cluster_assignments()[i][j] == _bg_cluster_tag) {
                                const seq_index_t index(i, j);
                                result += _precomputed_marginal[index];
                        }
                }
        }

        return result;
}

string
independence_background_t::print_counts() const {
        return string();
}

void
independence_background_t::set_bg_cluster_tag(cluster_tag_t bg_cluster_tag) {
        _bg_cluster_tag = bg_cluster_tag;
}

// Multinomial/Dirichlet Model
////////////////////////////////////////////////////////////////////////////////

product_dirichlet_t::product_dirichlet_t(
        const matrix<double>& _alpha,
        const sequence_data_t<data_tfbs_t::code_t>& data,
        const sequence_data_t<data_tfbs_t::code_t>& complement_data,
        bool complement)
        : component_model_t(),
          _size1(_alpha.size()),
          _size2(_alpha[0].size()),
          _data(&data),
          _complement_data(&complement_data),
          _complement(complement)
{
        for (size_t i = 0; i < _alpha.size(); i++) {
                alpha .push_back(counts_t());
                counts.push_back(counts_t());
                for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                        alpha [i][j] = _alpha[i][j];
                        counts[i][j] = _alpha[i][j];
                }
        }
}

product_dirichlet_t::product_dirichlet_t(const product_dirichlet_t& distribution)
        : component_model_t(distribution),
          alpha (distribution.alpha),
          counts(distribution.counts),
          _size1(distribution._size1),
          _size2(distribution._size2),
          _data (distribution._data),
          _complement_data(distribution._complement_data),
          _complement(distribution._complement)
{
}

product_dirichlet_t::~product_dirichlet_t() {
}

product_dirichlet_t*
product_dirichlet_t::clone() const {
        return new product_dirichlet_t(*this);
}

product_dirichlet_t&
product_dirichlet_t::operator=(const component_model_t& component_model)
{
        product_dirichlet_t tmp(static_cast<const product_dirichlet_t&>(component_model));
        swap(*this, tmp);
        return *this;
}

size_t
product_dirichlet_t::add(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();
        size_t i, k;

        if (!range.reverse() || _complement == false) {
                for (i = 0; i < length; i++) {
                        const seq_index_t index(sequence, position+i);
                        for (k = 0; k < data_tfbs_t::alphabet_size; k++) {
                                counts[i%_size1][k] += data()[index][k];
                        }
                }
        }
        // reverse complement
        else {
                for (i = 0; i < length; i++) {
                        const seq_index_t index(sequence, position+length-i-1);
                        for (k = 0; k < data_tfbs_t::alphabet_size; k++) {
                                counts[i%_size1][k] += complement_data()[index][k];
                        }
                }
        }
        return i/_size1;
}

size_t
product_dirichlet_t::remove(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();
        size_t i, k;

        if (!range.reverse() || _complement == false) {
                for (i = 0; i < length; i++) {
                        const seq_index_t index(sequence, position+i);
                        for (k = 0; k < data_tfbs_t::alphabet_size; k++) {
                                counts[i%_size1][k] -= data()[index][k];
                        }
                }
        }
        // reverse complement
        else {
                for (i = 0; i < length; i++) {
                        const seq_index_t index(sequence, position+length-i-1);
                        for (k = 0; k < data_tfbs_t::alphabet_size; k++) {
                                counts[i%_size1][k] -= complement_data()[index][k];
                        }
                }
        }
        return i/_size1;
}

size_t
product_dirichlet_t::count(const range_t& range) {
        return range.length()/_size1;
}

/*
 *  p(y|x) = Beta(n(x) + n(y) + alpha) / Beta(n(x) + alpha)
 */
double product_dirichlet_t::predictive(const range_t& range) {
        return exp(log_predictive(range));
}

double product_dirichlet_t::predictive(const vector<range_t>& range_set) {
        return exp(log_predictive(range_set));
}

double product_dirichlet_t::log_predictive(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();
        double result = 0;

        if (!range.reverse() || _complement == false) {
                for (size_t i = 0; i < length; i++) {
                        const seq_index_t index(sequence, position+i);

                        /* counts contains the data count statistic
                         * and the pseudo counts alpha */
                        result += fast_lnbeta<data_tfbs_t::alphabet_size>(counts[i%_size1], data()[index])
                                - fast_lnbeta<data_tfbs_t::alphabet_size>(counts[i%_size1]);
                }
        }
        // reverse complement
        else {
                for (size_t i = 0; i < length; i++) {
                        const seq_index_t index(sequence, position+length-i-1);

                        /* counts contains the data count statistic
                         * and the pseudo counts alpha */
                        result += fast_lnbeta<data_tfbs_t::alphabet_size>(counts[i%_size1], complement_data()[index])
                                - fast_lnbeta<data_tfbs_t::alphabet_size>(counts[i%_size1]);
                }
        }

        return result;
}

double product_dirichlet_t::log_predictive(const vector<range_t>& range_set) {
        assert(range_set.size() > 0);

        const size_t length = range_set[0].length();
        double result = 0;

        for (size_t i = 0; i < length; i++) {

                /* set all tmp_counts to zero */
                for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                        tmp_counts[j] = 0;
                }

                /* loop through all ranges */
                for (size_t k = 0; k < range_set.size(); k++) {

                        const size_t sequence = range_set[k].index()[0];
                        const size_t position = range_set[k].index()[1];
                        const seq_index_t index(sequence, position+i);

                        /* add counts of this subsequence to tmp_counts */
                        for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                                tmp_counts[j] += data()[index][j];
                        }
                }
                result += fast_lnbeta<data_tfbs_t::alphabet_size>(counts[i%_size1], tmp_counts)
                        - fast_lnbeta<data_tfbs_t::alphabet_size>(counts[i%_size1]);
        }

        return result;
}

/*
 *  p(x) = Beta(n(x) + alpha) / Beta(alpha)
 */
double product_dirichlet_t::log_likelihood() const {
        double result = 0;

        for (size_t i = 0; i < _size1; i++) {
                /* counts contains the data count statistic
                 * and the pseudo counts alpha */
                result += fast_lnbeta<data_tfbs_t::alphabet_size>(counts[i%_size1])
                        - fast_lnbeta<data_tfbs_t::alphabet_size>(alpha [i%_size1]);
        }
        return result;
}

string
product_dirichlet_t::print_counts() const {
        stringstream ss;
        for (size_t k = 0; k < _size2; k++) {
                ss << nucleotide_alphabet_t().decode(k) << " ";
                for (size_t j = 0; j < _size1; j++) {
                        ss.precision(8);
                        ss.width(10);
                        ss << counts[j][k] << " ";
                }
                ss << endl;
        }
        return ss.str();
}

// Markov Chain Mixture
////////////////////////////////////////////////////////////////////////////////

markov_chain_mixture_t::markov_chain_mixture_t(
        size_t alphabet_size,
        const tfbs_options_t& options,
        const sequence_data_t<data_tfbs_t::code_t>& data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        cluster_tag_t cluster_tag)
        : component_model_t(cluster_assignments),
          _data(data),
          _cluster_tag(cluster_tag),
          _max_context(options.background_context),
          _alphabet_size(alphabet_size)
{
//        cerr << "Error: The markov chain mixture model is broken!"
//             << endl;
//        exit(EXIT_FAILURE);

        _length = context_t::counts_size(alphabet_size, _max_context);
        _alpha  = (double*)malloc(_length*sizeof(double));
        _counts = (double*)malloc(_length*sizeof(double));
        _counts_sum = (double*)malloc(_length/_alphabet_size*sizeof(double));
        _parents = (int*)malloc(_length*sizeof(int));
        if (options.background_weights == "entropy") {
                _weights = new entropy_weights_t(_alphabet_size, _max_context, _length);
        }
        else if (options.background_weights == "decay") {
                _weights = new decay_weights_t(_max_context);
        }
        else {
                cerr << "Error: Unknown background weights." << endl;
                exit(EXIT_FAILURE);
        }

        /* for likelihood computations */
        _counts_tmp = (double*)malloc(_length*sizeof(double));

        /* init counts_sum to zero */
        for (size_t i = 0; i < _length/_alphabet_size; i++) {
                _counts_sum[i] = 0;
        }
        /* init counts */
        for (size_t i = 0; i < _length; i++) {
                _alpha[i]  = options.background_alpha[0][0];
                _counts[i] = 0.0;
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
                const sequence_t<double>& seq = (const sequence_t<double>&)_data[i];
                _context.push_back(seq_context_t(seq, _max_context, alphabet_size));
        }
}

markov_chain_mixture_t::markov_chain_mixture_t(const markov_chain_mixture_t& distribution)
        : component_model_t(distribution),
          _length(distribution._length),
          _weights(distribution._weights->clone()),
          _data(distribution._data),
          _context(distribution._context),
          _cluster_tag(distribution._cluster_tag),
          _max_context(distribution._max_context),
          _alphabet_size(distribution._alphabet_size)
{
        _alpha      = (double*)malloc(_length*sizeof(double));
        _counts     = (double*)malloc(_length*sizeof(double));
        _counts_sum = (double*)malloc(_length/_alphabet_size*sizeof(double));
        _parents    = (int*)malloc(_length*sizeof(int));

        /* for likelihood computations */
        _counts_tmp = (double*)malloc(_length*sizeof(double));

        /* init data structures */
        memcpy(_alpha,   distribution._alpha,   _length*sizeof(double));
        memcpy(_counts,  distribution._counts,  _length*sizeof(double));
        memcpy(_counts_sum, distribution._counts_sum, _length/_alphabet_size*sizeof(double));
        memcpy(_parents, distribution._parents, _length*sizeof(int));
}

markov_chain_mixture_t::~markov_chain_mixture_t() {
        free(_alpha);
        free(_counts);
        free(_parents);
        free(_counts_sum);
        free(_counts_tmp);
        delete(_weights);
}

markov_chain_mixture_t*
markov_chain_mixture_t::clone() const {
        return new markov_chain_mixture_t(*this);
}

markov_chain_mixture_t&
markov_chain_mixture_t::operator=(const component_model_t& component_model)
{
        // broken
        exit(EXIT_FAILURE);
        return *this;
}

size_t
markov_chain_mixture_t::max_from_context(const range_t& range) const
{
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        ssize_t from = position;

        for (size_t i = 0; i < _max_context && from > 0 && cluster_assignments()[seq_index_t(sequence, from-1)] == _cluster_tag; i++) {
                from--;
        }

        return position-from;
}

size_t
markov_chain_mixture_t::max_to_context(const range_t& range) const
{
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();
        const ssize_t sequence_length = _data.size(sequence);
        ssize_t to = position+length-1;

        for (size_t i = 0; i < _max_context && to < sequence_length-1 && cluster_assignments()[seq_index_t(sequence, to + 1)] == _cluster_tag; i++) {
                to++;
        }

        return to-position-length+1;
}

size_t
markov_chain_mixture_t::add(const range_t& range) {
        const size_t sequence     = range.index()[0];
        const size_t length       = range.length();
        const size_t from_context = max_from_context(range);
        const size_t   to_context = max_to_context(range);

        for (size_t i = 0; i < length+to_context; i++) {
                const size_t pos   = range.index()[1]+i;
                const size_t c_max = min(from_context+i, _max_context);
                const size_t c_min = i < length ? 0 : i-(length-1);
                for (size_t c = c_min; c <= c_max; c++) {
                        const int code = _context[sequence][pos][c];
                        if (code != -1) {
                                _counts[code]++;
                                _counts_sum[code/_alphabet_size]++;
                                _weights->update(code, _alpha, _counts, _counts_sum);
                        }
                }
        }

        return range.length();
}

size_t
markov_chain_mixture_t::remove(const range_t& range) {
        const size_t sequence     = range.index()[0];
        const size_t length       = range.length();
        const size_t from_context = max_from_context(range);
        const size_t   to_context = max_to_context(range);

        for (size_t i = 0; i < length+to_context; i++) {
                const size_t pos   = range.index()[1]+i;
                const size_t c_max = min(from_context+i, _max_context);
                const size_t c_min = i < length ? 0 : i-(length-1);
                for (size_t c = c_min; c <= c_max; c++) {
                        const int code = _context[sequence][pos][c];
                        if (code != -1) {
                                _counts[code]--;
                                _counts_sum[code/_alphabet_size]--;
                                _weights->update(code, _alpha, _counts, _counts_sum);
                        }
                }
        }

        return range.length();
}

size_t
markov_chain_mixture_t::count(const range_t& range) {
        return range.length();
}

double markov_chain_mixture_t::predictive(const range_t& range) {
        const size_t sequence     = range.index()[0];
        const size_t length       = range.length();
        const size_t from_context = max_from_context(range);
        const size_t   to_context = max_to_context(range);
        double result = 1;

        for (size_t i = 0; i < length+to_context; i++) {
                const size_t pos   = range.index()[1]+i;
                const size_t c_max = min(from_context+i, _max_context);
                const size_t c_min = i < length ? 0 : i-(length-1);
                vector<int> codes(c_max-c_min+1, 0);
                for (size_t c = c_min; c <= c_max; c++) {
                        codes[c-c_min] = _context[sequence][pos][c];
                }
                double partial_result = 0;
                _weights->init(codes);
                for (size_t c = c_min; c <= c_max; c++) {
                        const int code = _context[sequence][pos][c];
                        // compute mixture component
                        partial_result +=
                                (*_weights)[c-c_min]*
                                (_counts[code]+_alpha[code])/
                                _counts_sum[code/_alphabet_size];
                }
                // save result
                result *= partial_result;
        }

        return result;
}

double markov_chain_mixture_t::predictive(const vector<range_t>& range) {
        assert(false);
}

double markov_chain_mixture_t::log_predictive(const range_t& range) {
        return log(predictive(range));
}

double markov_chain_mixture_t::log_predictive(const vector<range_t>& range_set) {
        return log(predictive(range_set));
}

double markov_chain_mixture_t::log_likelihood(size_t pos) const
{
        const double n = _counts_tmp[pos];
        double result = 0;
        vector<int> path;
        vector<int> path_rev;

        // compute the path up the tree
        for (int i = pos; i != -1; i = _parents[i]) {
                path.push_back(i);
        }
        reverse(path.begin(), path.end());
        // compute likelihood for this node
        _weights->init(path);
        for (size_t c = 0; c < path.size(); c++) {
                const int code = path[c];
                // compute mixture component
                result += (*_weights)[c]*
                        (_counts[code]+_alpha[code])/
                        _counts_sum[code/_alphabet_size];
        }

        return n*log(result);
}

void markov_chain_mixture_t::substract_counts(size_t pos) const
{
        const double n = _counts_tmp[pos];

        for (int i = pos; i != -1; i = _parents[i]) {
                _counts_tmp[i] -= n;
        }
}

double markov_chain_mixture_t::log_likelihood() const {
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

// Bivariate Gaussian
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
