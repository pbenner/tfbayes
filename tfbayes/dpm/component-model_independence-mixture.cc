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

#include <algorithm>    // std::max_element
#include <cmath>        // std::exp, std::log
#include <string>
#include <vector>

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/fastarithmetics/fast-lnbeta.hh>

using namespace std;

// Independence Model
////////////////////////////////////////////////////////////////////////////////

/* This is an independence background model with Dirichlet
 * prior. Each position in the alignment is drawn independently from
 * a mixture of n components, where n is determined by
 * the number of columns of alpha. Alignment columns are assigned to
 * a single component by maximizing the posterior probability. */
independence_mixture_background_t::independence_mixture_background_t(
        const matrix<double>& _alpha,
        const vector<double>& weights,
        const sequence_data_t<data_tfbs_t::code_t>& _data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        boost::optional<const alignment_set_t<>&> alignment_set)
        : component_model_t({"background", 1}, cluster_assignments),
          _size(data_tfbs_t::alphabet_size),
          _bg_cluster_tag(0),
          _precomputed_marginal(_data.sizes(), 0),
          _data(&_data)
{
        vector<counts_t> alpha(_alpha.size(), counts_t());
        // check the pseudocounts matrix (each column is a component
        // of the background model)
        assert(_alpha.size() > 0);
        assert(_alpha.size() == weights.size());
        for (size_t i = 0; i < alpha.size(); i++) {
                assert(_alpha[i].size() == data_tfbs_t::alphabet_size);
                for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                        alpha[i][j] = _alpha[i][j];
                }
        }
        precompute_marginal(alpha, weights);
}

independence_mixture_background_t::independence_mixture_background_t(const independence_mixture_background_t& distribution)
        : component_model_t(distribution),
          _size(distribution._size),
          _bg_cluster_tag(distribution._bg_cluster_tag),
          _precomputed_marginal(distribution._precomputed_marginal),
          _data(distribution._data)
{
}

independence_mixture_background_t::~independence_mixture_background_t() {
}

independence_mixture_background_t*
independence_mixture_background_t::clone() const {
        return new independence_mixture_background_t(*this);
}

independence_mixture_background_t&
independence_mixture_background_t::operator=(const component_model_t& component_model)
{
        independence_mixture_background_t tmp(
                static_cast<const independence_mixture_background_t&>(component_model));
        swap(*this, tmp);
        return *this;
}

void
independence_mixture_background_t::precompute_marginal(
        const vector<counts_t>& alpha,
        const vector<double  >& weights)
{
        // marginals for each component
        vector<double> tmp(alpha.size(), 0.0);
        // keep track of the number of positions assigned to each
        // component
        vector<size_t> statistics(alpha.size(), 0);
        size_t n = 0;
        /* go through the data and precompute
         * lnbeta(n + alpha) - lnbeta(alpha) */
        for(size_t i = 0; i < data().size(); i++) {
                for(size_t j = 0; j < data()[i].size(); j++) {
                        for (size_t c = 0; c < alpha.size(); c++) {
                                tmp[c] = fast_lnbeta(alpha[c], data()[i][j])
                                       - fast_lnbeta(alpha[c])
                                       + log(weights[c]);
                        }
                        size_t c = distance(tmp.begin(),
                                            max_element(tmp.begin(), tmp.end()));
                        _precomputed_marginal[i][j] = tmp[c] - log(weights[c]);
                        // increment count statistics
                        statistics[c]++; n++;
                }
        }
        flockfile(stderr);
        cerr << "Background assignments:"
             << endl;
        for (size_t c = 0; c < alpha.size(); c++) {
                cerr << boost::format("-> %3.2f%% assigned to component %d\n")
                        % (100.0*statistics[c]/n) % (c+1);
        }
        funlockfile(stderr);        
}

size_t
independence_mixture_background_t::add(const range_t& range) {
        return range.length();
}

size_t
independence_mixture_background_t::remove(const range_t& range) {
        return range.length();
}

size_t
independence_mixture_background_t::count(const range_t& range) {
        return range.length();
}

/*
 *  p(y|x) = Beta(n(x) + n(y) + alpha) / Beta(n(x) + alpha)
 */
double independence_mixture_background_t::predictive(const range_t& range) {
        return exp(log_predictive(range));
}

double independence_mixture_background_t::predictive(const vector<range_t>& range_set) {
        return exp(log_predictive(range_set));
}

double independence_mixture_background_t::log_predictive(const range_t& range) {
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

double independence_mixture_background_t::log_predictive(const vector<range_t>& range_set) {
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
double independence_mixture_background_t::log_likelihood() const {
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
independence_mixture_background_t::print_counts() const {
        return string();
}

void
independence_mixture_background_t::set_bg_cluster_tag(cluster_tag_t bg_cluster_tag) {
        _bg_cluster_tag = bg_cluster_tag;
}
