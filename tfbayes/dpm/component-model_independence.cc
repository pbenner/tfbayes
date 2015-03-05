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

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

#include <boost/format.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/fastarithmetics/fast-lnbeta.hh>

using namespace std;

// Methods to integrate out the pseudocounts of a Dirichlet
// distribution, where we use a Gamma hyperprior
////////////////////////////////////////////////////////////////////////////////

struct gamma_marginal_data {
        const independence_background_t::counts_t& counts;
        const independence_background_t::counts_t& alpha;
        const boost::math::gamma_distribution<> distribution;
};

GCC_ATTRIBUTE_HOT
double
gamma_marginal_f(double * x, size_t dim, void * params)
{
        /* store the result on normal scale */
        double result = 1.0;
        /* casted parameters */
        struct gamma_marginal_data* data = (gamma_marginal_data *)params;

        /* pseudocounts */
        independence_background_t::counts_t alpha(data->alpha);

        for (size_t i = 0, j = 0; i < data_tfbs_t::alphabet_size; i++) {
                /* determine pseudocounts that are integrated out */
                if (alpha[i] == -1) {
                        assert(j < dim);
                        /* copy pseudocounts */
                        alpha[i] = x[j];
                        /* multiply with gamma distribution */
                        result  *= boost::math::pdf(data->distribution, x[j++]);
                }
        }

        /* lnbeta(n, a) - lnbeta(a) */
        result *= exp(fast_lnbeta(data->counts, alpha) -
                      fast_lnbeta(alpha));

        return result;
}

/* use a gamma prior distribution for the Dirichlet pseudocounts to
 * integrate them out (there is no closed form solution so we need to
 * do this numerically) */
double
gamma_marginal(
        const independence_background_t::counts_t& counts,
        const independence_background_t::counts_t& alpha,
        const vector<double>& parameters)
{
        size_t dim = count(alpha.begin(), alpha.end(), -1);
        double xl[data_tfbs_t::alphabet_size];
        double xu[data_tfbs_t::alphabet_size];
        double k = parameters[0];
        double g = parameters[1];
        const gsl_rng_type *T;
        gsl_rng *r;

        size_t calls = 500000;
        double result, err;

        gsl_monte_function F;

        struct gamma_marginal_data data = {
                counts, alpha, boost::math::gamma_distribution<>(
                        parameters[0], parameters[1])
        };

        /* begin at the mode and determine a point where the density
         * function is below a certain threshold */
        double thr;
        for (thr = (k-1)*g; boost::math::pdf(data.distribution, thr) > 1e-8; thr += 1.0);

        for (size_t i = 0; i < dim; i++) {
                xl[i] = 0.0;
                xu[i] = thr;
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
        const vector<double>& parameters,
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
                double result = gamma_marginal(counts, alpha, parameters);
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
                ar & parameters;
                ar & data;
                ar & precomputed_marginal;
        }
public:
        background_cache_t()
                { }
        background_cache_t(
                const independence_background_t::counts_t& alpha,
                const vector<double>& parameters,
                const sequence_data_t<data_tfbs_t::code_t>& data,
                const sequence_data_t<double>& precomputed_marginal)
                : alpha(alpha), parameters(parameters), data(data),
                  precomputed_marginal(precomputed_marginal)
                { }
        bool consistent(
                const independence_background_t::counts_t& alpha,
                const vector<double>& parameters,
                const sequence_data_t<data_tfbs_t::code_t>& data) {
                // check pseudocounts
                if (this->alpha.size() != alpha.size()) {
                        return false;
                }
                for (size_t i = 0; i < alpha.size(); i++) {
                        if (this->alpha[i] != alpha[i]) {
                                return false;
                        }
                }
                // gamma distribution parameters
                if (this->parameters.size() != parameters.size()) {
                        return false;
                }
                for (size_t i = 0; i < parameters.size(); i++) {
                        if (this->parameters[i] != parameters[i]) {
                                return false;
                        }
                }
                // check phylogenetic data
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
        vector<double> parameters;
        sequence_data_t<data_tfbs_t::code_t> data;
        sequence_data_t<double> precomputed_marginal;
};

// Independence Model
////////////////////////////////////////////////////////////////////////////////

/* This is an independence background model with Dirichlet
 * prior. If the pseudocounts are set to -1 a Gamma distribution
 * is used to integrate them out. */
independence_background_t::independence_background_t(
        const vector<double>& _alpha,
        const vector<double>& parameters,
        const sequence_data_t<data_tfbs_t::code_t>& _data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        thread_pool_t& thread_pool,
        const string& cachefile,
        boost::optional<const alignment_set_t<>&> alignment_set)
        : component_model_t({"background", 1}, cluster_assignments),
          _size(data_tfbs_t::alphabet_size),
          _bg_cluster_tag(0),
          _precomputed_marginal(_data.sizes(), 0),
          _data(&_data)
{
        assert(_alpha.size() == data_tfbs_t::alphabet_size);

        counts_t alpha;

        for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                alpha[j] = _alpha[j];
        }

        if (std::find(alpha.begin(), alpha.end(), -1) == alpha.end()) {
                // all pseudocounts are given
                precompute_marginal(alpha);
        }
        else {
                // at least one pseudocount has to be integrated
                if (!load_marginal_gamma(alpha, parameters, cachefile)) {
                        precompute_marginal_gamma(alpha, parameters, thread_pool);
                        save_marginal_gamma(alpha, parameters, cachefile);
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
        const vector<double>& parameters,
        const string& cachefile)
{
        std::ifstream ifs(cachefile, std::ios::binary);
        background_cache_t background_cache;

        // try to receive precomputed marginals from cache
        if (ifs) {
                boost::archive::binary_iarchive ia(ifs);
                ia >> background_cache;
                if (background_cache.consistent(alpha, parameters, data())) {
                        _precomputed_marginal = background_cache.precomputed_marginal;
                        return true;
                }
                flockfile(stderr);
                cerr << "Background cache is inconsistent... recomputing!"
                     << endl;
                funlockfile(stderr);
        }
        return false;
}

bool
independence_background_t::save_marginal_gamma(
        const counts_t& alpha,
        const vector<double>& parameters,
        const string& cachefile)
{
        std::ofstream ofs(cachefile);
        if (ofs) {
                boost::archive::binary_oarchive oa(ofs);
                const background_cache_t tmp(alpha, parameters, data(), _precomputed_marginal);
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
                                  fast_lnbeta(alpha, data()[i][j])
                                - fast_lnbeta(alpha);
                }
        }
}

void
independence_background_t::precompute_marginal_gamma(
        const counts_t& alpha,
        const vector<double>& parameters,
        thread_pool_t& thread_pool)
{
        typedef double (*hgm)(const counts_t&, const counts_t&,
                              const vector<double>&,
                              boost::unordered_map<counts_t, double>&,
                              boost::shared_mutex&);
        boost::unordered_map<counts_t, double> map;
        boost::shared_mutex mutex;

        flockfile(stderr);
        cerr << "Background gamma shape: " << parameters[0] << endl
             << "Background gamma scale: " << parameters[1] << endl;
        funlockfile(stderr);
        /* go through the data and precompute
         * lnbeta(n + alpha) - lnbeta(alpha) */
        for(size_t i = 0; i < data().size(); i++) {
                future_vector_t<double> futures(data()[i].size());

                for(size_t j = 0; j < data()[i].size(); j++) {
                        boost::function<double ()> f = boost::bind(
                                static_cast<hgm>(&hashed_gamma_marginal),
                                boost::cref(data()[i][j]),
                                boost::cref(alpha),
                                boost::cref(parameters),
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
        cerr << "\rPrecomputing background...   done." << endl << flush;
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
                const index_t index(sequence, position+i);

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
                        const index_t index(sequence, position+i);

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
                                const index_t index(i, j);
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
