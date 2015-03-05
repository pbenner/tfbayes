/* Copyright (C) 2011-2015 Philipp Benner
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

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/fastarithmetics/fast-lnbeta.hh>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

product_dirichlet_t::product_dirichlet_t(
        const model_id_t& model_id,
        const matrix<double>& _alpha,
        const sequence_data_t<data_tfbs_t::code_t>& data,
        const sequence_data_t<data_tfbs_t::code_t>& complement_data)
        : component_model_t(model_id),
          _size1(model_id.length),
          _size2(data_tfbs_t::alphabet_size),
          _data(&data),
          _complement_data(&complement_data)
{
        // make sure the counts vector has a single column
        assert(_alpha.size() == 1);
        assert(_alpha[0].size() == data_tfbs_t::alphabet_size);

        for (size_t i = 0; i < _size1; i++) {
                alpha .push_back(counts_t());
                counts.push_back(counts_t());
                for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                        alpha [i][j] = _alpha[0][j];
                        counts[i][j] = _alpha[0][j];
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
          _complement_data(distribution._complement_data)
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

        assert(length == _size1);

        if (!range.reverse()) {
                for (size_t i = 0; i < length; i++) {
                        const index_t index(sequence, position+i);
                        for (size_t k = 0; k < data_tfbs_t::alphabet_size; k++) {
                                counts[i][k] += data()[index][k];
                        }
                }
        }
        // reverse complement
        else {
                for (size_t i = 0; i < length; i++) {
                        const index_t index(sequence, position+length-i-1);
                        for (size_t k = 0; k < data_tfbs_t::alphabet_size; k++) {
                                counts[i][k] += complement_data()[index][k];
                        }
                }
        }
        return 1;
}

size_t
product_dirichlet_t::remove(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();

        assert(length == _size1);

        if (!range.reverse()) {
                for (size_t i = 0; i < length; i++) {
                        const index_t index(sequence, position+i);
                        for (size_t k = 0; k < data_tfbs_t::alphabet_size; k++) {
                                counts[i][k] -= data()[index][k];
                                assert(counts[i][k] >= 0.0);
                        }
                }
        }
        // reverse complement
        else {
                for (size_t i = 0; i < length; i++) {
                        const index_t index(sequence, position+length-i-1);
                        for (size_t k = 0; k < data_tfbs_t::alphabet_size; k++) {
                                counts[i][k] -= complement_data()[index][k];
                                assert(counts[i][k] >= 0.0);
                        }
                }
        }
        return 1;
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

        // return zero if the length does not match
        if (length != _size1) {
                return -std::numeric_limits<double>::infinity();
        }
        if (!range.reverse()) {
                for (size_t i = 0; i < length; i++) {
                        const index_t index(sequence, position+i);

                        /* counts contains the data count statistic
                         * and the pseudo counts alpha */
                        result += fast_lnbeta(counts[i], data()[index])
                                - fast_lnbeta(counts[i]);
                }
        }
        // reverse complement
        else {
                for (size_t i = 0; i < length; i++) {
                        const index_t index(sequence, position+length-i-1);

                        /* counts contains the data count statistic
                         * and the pseudo counts alpha */
                        result += fast_lnbeta(counts[i], complement_data()[index])
                                - fast_lnbeta(counts[i]);
                }
        }

        return result;
}

double product_dirichlet_t::log_predictive(const vector<range_t>& range_set) {
        assert(range_set.size() > 0);

        const size_t length = range_set[0].length();
        double result = 0;

        assert(length == _size1);

        for (size_t i = 0; i < length; i++) {

                /* set all tmp_counts to zero */
                for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                        tmp_counts[j] = 0;
                }

                /* loop through all ranges */
                for (size_t k = 0; k < range_set.size(); k++) {

                        if (!range_set[k].reverse()) {
                                const size_t sequence = range_set[k].index()[0];
                                const size_t position = range_set[k].index()[1];
                                const index_t index(sequence, position+i);

                                /* add counts of this subsequence to tmp_counts */
                                for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                                        tmp_counts[j] += data()[index][j];
                                }
                        }
                }
                /* reverse complements */
                for (size_t k = 0; k < range_set.size(); k++) {

                        if (range_set[k].reverse()) {
                                const size_t sequence = range_set[k].index()[0];
                                const size_t position = range_set[k].index()[1];
                                const index_t index(sequence, position+length-i-1);

                                /* add counts of this subsequence to tmp_counts */
                                for (size_t j = 0; j < data_tfbs_t::alphabet_size; j++) {
                                        tmp_counts[j] += complement_data()[index][j];
                                }
                        }
                }
                result += fast_lnbeta(counts[i], tmp_counts)
                        - fast_lnbeta(counts[i]);
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
                result += fast_lnbeta(counts[i])
                        - fast_lnbeta(alpha [i]);
        }
        return result;
}

string
product_dirichlet_t::print_counts() const {
        stringstream ss;
        for (size_t k = 0; k < _size2; k++) {
                ss << nucleotide_alphabet_t().decode(k) << " ";
                for (size_t j = 0; j < _size1; j++) {
                        ss.width(10);
                        ss.precision(2);
                        ss << fixed << counts[j][k] << " ";
                }
                ss << endl;
        }
        return ss.str();
}
