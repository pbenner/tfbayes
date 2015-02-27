/* Copyright (C) 2015 Philipp Benner
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

#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/utility/probability.hh>
#include <tfbayes/entropy/entropy-multinomial-distribution.hh>
#include <tfbayes/fastarithmetics/fast-lnbeta.hh>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

class background_cache_t {
private:
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
                ar & id;
                ar & parameters;
                ar & data;
                ar & precomputed_marginal;
        }
public:
        background_cache_t()
                : id("entropy")
                { }
        background_cache_t(
                const vector<double>& parameters,
                const sequence_data_t<data_tfbs_t::code_t>& data,
                const sequence_data_t<double>& precomputed_marginal)
                : id         ("entropy")
                , parameters (parameters)
                , data       (data)
                , precomputed_marginal(precomputed_marginal)
                { }
        bool consistent(
                const vector<double>& parameters,
                const sequence_data_t<data_tfbs_t::code_t>& data) {
                // id
                if (this->id != "entropy") {
                        return false;
                }
                // parameters
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

        string id;
        vector<double> parameters;
        sequence_data_t<data_tfbs_t::code_t> data;
        sequence_data_t<double> precomputed_marginal;
};

////////////////////////////////////////////////////////////////////////////////

/* This is an independence background model with Dirichlet
 * prior. If the pseudocounts are set to -1 a Gamma distribution
 * is used to integrate them out. */
entropy_background_t::entropy_background_t(
        const vector<double>& parameters,
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
        assert(parameters.size() == 2);

        // at least one pseudocount has to be integrated
        if (!load_marginal(parameters, cachefile)) {
                precompute_marginal(parameters, thread_pool);
                save_marginal(parameters, cachefile);
        }
}

entropy_background_t::entropy_background_t(const entropy_background_t& distribution)
        : component_model_t(distribution),
          _size(distribution._size),
          _bg_cluster_tag(distribution._bg_cluster_tag),
          _precomputed_marginal(distribution._precomputed_marginal),
          _data(distribution._data)
{
}

entropy_background_t::~entropy_background_t() {
}

entropy_background_t*
entropy_background_t::clone() const {
        return new entropy_background_t(*this);
}

entropy_background_t&
entropy_background_t::operator=(const component_model_t& component_model)
{
        entropy_background_t tmp(
                static_cast<const entropy_background_t&>(component_model));
        swap(*this, tmp);
        return *this;
}

bool
entropy_background_t::load_marginal(
        const vector<double>& parameters,
        const string& cachefile)
{
        std::ifstream ifs(cachefile, std::ios::binary);
        background_cache_t background_cache;

        // try to receive precomputed marginals from cache
        if (ifs) {
                boost::archive::binary_iarchive ia(ifs);
                ia >> background_cache;
                if (background_cache.consistent(parameters, data())) {
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
entropy_background_t::save_marginal(
        const vector<double>& parameters,
        const string& cachefile)
{
        std::ofstream ofs(cachefile);
        if (ofs) {
                boost::archive::binary_oarchive oa(ofs);
                const background_cache_t tmp(parameters, data(), _precomputed_marginal);
                oa << tmp;
                flockfile(stderr);
                cerr << boost::format("Background cache saved to `%s'.") % cachefile
                     << endl;
                funlockfile(stderr);

                return true;
        }
        return false;
}

class precompute_marginal_functor
{
        template <typename T>
        void seed_rng(T& rng) {
                struct timeval tv;
                gettimeofday(&tv, NULL);
                rng.seed(tv.tv_sec*tv.tv_usec);
        }
        typedef double real_t;
        // the marginal entropy distribution requires a special
        // probability type
        typedef probability_t<real_t> p_t;
public:
        typedef real_t result_type;

        precompute_marginal_functor(const vector<real_t>& parameters) {

                boost::random::mt19937 gen;
                seed_rng(gen);

                m_dist = marginal_entropy_distribution_t<double, p_t>(
                        data_tfbs_t::alphabet_size, parameters[0], parameters[1], 100000, gen);
        }

        template <class county_type>
        result_type operator()(const county_type& counts) const {
                return log(pdf(m_dist, counts));
        }
protected:
        marginal_entropy_distribution_t<real_t, p_t> m_dist;
};

void
entropy_background_t::precompute_marginal(
        const vector<double>& parameters,
        thread_pool_t& thread_pool)
{
        precompute_marginal_functor functor(parameters);

        flockfile(stderr);
        cerr << boost::format("Background beta pseudocounts: %f, %f")
                              % parameters[0] % parameters[1]
             << endl;
        funlockfile(stderr);

        // go through the data and precompute the marginal distribution
        for(size_t i = 0; i < data().size(); i++) {
                future_vector_t<double> futures(data()[i].size());

                for(size_t j = 0; j < data()[i].size(); j++) {
                        boost::function<double ()> f = boost::bind(
                                boost::cref(functor), boost::cref(data()[i][j]));

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
entropy_background_t::add(const range_t& range) {
        return range.length();
}

size_t
entropy_background_t::remove(const range_t& range) {
        return range.length();
}

size_t
entropy_background_t::count(const range_t& range) {
        return range.length();
}

/*
 *  p(y|x) = Beta(n(x) + n(y) + alpha) / Beta(n(x) + alpha)
 */
double entropy_background_t::predictive(const range_t& range) {
        return exp(log_predictive(range));
}

double entropy_background_t::predictive(const vector<range_t>& range_set) {
        return exp(log_predictive(range_set));
}

double entropy_background_t::log_predictive(const range_t& range) {
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

double entropy_background_t::log_predictive(const vector<range_t>& range_set) {
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
double entropy_background_t::log_likelihood() const {
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
entropy_background_t::print_counts() const {
        return string();
}

void
entropy_background_t::set_bg_cluster_tag(cluster_tag_t bg_cluster_tag) {
        _bg_cluster_tag = bg_cluster_tag;
}
