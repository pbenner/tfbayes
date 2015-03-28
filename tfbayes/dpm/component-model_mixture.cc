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

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/utility/normalize.hh>
#include <tfbayes/utility/multinomial-beta.hh>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
mixture_dirichlet_t::mixture_dirichlet_t(
        const matrix<double>& _alpha,
        const vector<double>& _weights,
        const sequence_data_t<data_tfbs_t::code_t>& data)
        : component_model_t({"background", 1}),
          weights(normalize(_weights)),
          _size1 (_alpha.size()),
          _size2 (_alpha[0].size()),
          _data  (&data),
          _component_assignments(data.sizes(), -1)
{
        if (weights.size() == 0) {
                weights = vector<double>(_size1, 1.0/static_cast<double>(_size1));
        }
        else {
                assert(weights.size() == _size1);
        }
        for (size_t i = 0; i < _alpha.size(); i++) {
                alpha .push_back(counts_t());
                counts.push_back(counts_t());
                for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                        alpha [i][j] = _alpha[i][j];
                        counts[i][j] = _alpha[i][j];
                }
        }
        precompute_component_assignments();
}

ssize_t
mixture_dirichlet_t::max_component(const index_t& index) const
{
        vector<double> result(_size1, 0.0);

        /* find component with highest predictive value */
        for (size_t i = 0; i < _size1; i++) {
                /* counts contains the data count statistic
                 * and the pseudo counts alpha */
                result[i] = mbeta_log(counts[i], data()[index])
                          - mbeta_log(counts[i])
                          + log(weights[i]);
        }
        return distance(result.begin(), max_element(result.begin(), result.end()));
}

bool
mixture_dirichlet_t::precompute_component_assignments_loop()
{
        bool optimized = false;
        /* optimize assignments */
        for (size_t i = 0; i < data().size(); i++) {
                for (size_t j = 0; j < data()[i].size(); j++) {
                        index_t index(i, j);
                        /* remove counts if possible */
                        if (_component_assignments[index] != -1) {
                                remove(index);
                        }
                        /* get best assignment */
                        ssize_t k = max_component(index);
                        /* check if assigment changed */
                        if (k != _component_assignments[index]) {
                                optimized = true;
                        }
                        /* save assignemnt */
                        _component_assignments[index] = k;
                        /* add counts to background */
                        add(index);
                }
        }
        return optimized;
}

#include <tfbayes/utility/terminal-codes.hh>

void
mixture_dirichlet_t::precompute_component_assignments()
{
        bool optimized = true;
        /* iterate until reaching a fixed point */
        flockfile(stderr);
        cerr << endl;
        for (size_t i = 0; optimized; i++) {
                cerr << __line_up__ << __line_del__
                     << boost::format("Optimizing background assignments ... [%d]") % i
                     << endl;

                optimized = precompute_component_assignments_loop();
        }
        cerr << __line_up__ << __line_del__
             << "Optimizing background assignments ... done"
             << endl;
        funlockfile(stderr);
}

mixture_dirichlet_t::mixture_dirichlet_t(const mixture_dirichlet_t& distribution)
        : component_model_t(distribution),
          alpha   (distribution.alpha),
          counts  (distribution.counts),
          weights (distribution.weights),
          _size1  (distribution._size1),
          _size2  (distribution._size2),
          _data   (distribution._data),
          _component_assignments(distribution._component_assignments)
{
}

mixture_dirichlet_t::~mixture_dirichlet_t() {
}

mixture_dirichlet_t*
mixture_dirichlet_t::clone() const {
        return new mixture_dirichlet_t(*this);
}

mixture_dirichlet_t&
mixture_dirichlet_t::operator=(const component_model_t& component_model)
{
        mixture_dirichlet_t tmp(static_cast<const mixture_dirichlet_t&>(component_model));
        swap(*this, tmp);
        return *this;
}

size_t
mixture_dirichlet_t::add(const index_t& index)
{
        size_t i = _component_assignments[index];

        /* add counts to the max component */
        for (size_t k = 0; k < data_tfbs_t::alphabet_size; k++) {
                counts[i][k] += data()[index][k];
        }

        return 1;
}

size_t
mixture_dirichlet_t::add(const range_t& range)
{
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();

        for (size_t i = 0; i < length; i++) {
                add(index_t(sequence, position+i));
        }
        return length;
}

size_t
mixture_dirichlet_t::remove(const index_t& index)
{
        size_t i = _component_assignments[index];

        /* substract counts from the max component */
        for (size_t k = 0; k < data_tfbs_t::alphabet_size; k++) {
                counts[i][k] -= data()[index][k];
        }

        return 1;
}

size_t
mixture_dirichlet_t::remove(const range_t& range)
{
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();

        for (size_t i = 0; i < length; i++) {
                remove(index_t(sequence, position+i));
        }
        return length;
}

size_t
mixture_dirichlet_t::count(const range_t& range) {
        return range.length();
}

/*
 *  p(y|x) = Beta(n(x) + n(y) + alpha) / Beta(n(x) + alpha)
 */
double mixture_dirichlet_t::predictive(const range_t& range) {
        return exp(log_predictive(range));
}

double mixture_dirichlet_t::predictive(const vector<range_t>& range_set) {
        return exp(log_predictive(range_set));
}

double mixture_dirichlet_t::log_predictive(const index_t& index) {
        size_t i = _component_assignments[index];

        /* counts contains the data count statistic
         * and the pseudo counts alpha */
        double result = mbeta_log(counts[i], data()[index])
                      - mbeta_log(counts[i]);

        return result;
}

double mixture_dirichlet_t::log_predictive(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();
        double result = 0;

        for (size_t i = 0; i < length; i++) {
                result += log_predictive(index_t(sequence, position+i));
        }
        // use high precision in debug mode
#ifdef DEBUG
        result = 0;
        add(range);
        result += log_likelihood();
        remove(range);
        result -= log_likelihood();
#endif
        return result;
}

double mixture_dirichlet_t::log_predictive(const vector<range_t>& range_set) {
        assert(range_set.size() > 0);
        double result = 0;

        BOOST_FOREACH(const range_t& range, range_set) {
                result += log_predictive(range);
        }
        // use high precision in debug mode
#ifdef DEBUG
        result = 0;
        BOOST_FOREACH(const range_t& range, range_set) {
                add(range);
        }
        result += log_likelihood();
        BOOST_FOREACH(const range_t& range, range_set) {
                remove(range);
        }
        result -= log_likelihood();
#endif

        return result;
}

/*
 *  p(x) = Beta(n(x) + alpha) / Beta(alpha)
 */
double mixture_dirichlet_t::log_likelihood() const {
        double result = 0;

        for (size_t i = 0; i < _size1; i++) {
                /* counts contains the data count statistic
                 * and the pseudo counts alpha */
                result += mbeta_log(counts[i])
                        - mbeta_log(alpha [i]);
        }
        return result;
}

string
mixture_dirichlet_t::print_counts() const {
        stringstream ss;
        for (size_t k = 0; k < _size2; k++) {
                ss << " -> " << nucleotide_alphabet_t().decode(k) << " ";
                for (size_t j = 0; j < _size1; j++) {
                        ss.width(10);
                        ss.precision(2);
                        ss << fixed << counts[j][k] << " ";
                }
                if (k+1 != _size2) {
                        ss << endl;
                }
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
        : component_model_t({"background", 1}, cluster_assignments),
          _data(data),
          _cluster_tag(cluster_tag),
          _max_context(options.background_context),
          _alphabet_size(alphabet_size)
{
       cerr << "Error: The markov chain mixture model is broken!"
            << endl;
       exit(EXIT_FAILURE);

        _length = context_t::counts_size(alphabet_size, _max_context);
        _alpha  = (double*)malloc(_length*sizeof(double));
        _counts = (double*)malloc(_length*sizeof(double));
        _counts_sum = (double*)malloc(_length/_alphabet_size*sizeof(double));
        _parents = (int*)malloc(_length*sizeof(int));
        _weights = new entropy_weights_t(_alphabet_size, _max_context, _length);
        //_weights = new decay_weights_t(_max_context);

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

        for (size_t i = 0; i < _max_context && from > 0 && cluster_assignments()[index_t(sequence, from-1)] == _cluster_tag; i++) {
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

        for (size_t i = 0; i < _max_context && to < sequence_length-1 && cluster_assignments()[index_t(sequence, to + 1)] == _cluster_tag; i++) {
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
