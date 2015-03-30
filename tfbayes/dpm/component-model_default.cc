/* Copyright (C) 2014 Philipp Benner
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
#include <tfbayes/utility/multinomial-beta.hh>

#include <boost/math/distributions/gamma_extra.hpp>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

default_background_t::default_background_t(
        const vector<double>& parameters,
        const sequence_data_t<data_tfbs_t::code_t>& _data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        thread_pool_t& thread_pool,
        const string& cachefile,
        boost::optional<const alignment_set_t<>&> alignment_set,
        size_t verbose)
        : component_model_t      ({"background", 1}, cluster_assignments)
        , prior_distribution     (parameters[0], parameters[1])
        , m_size                 (data_tfbs_t::alphabet_size)
        , m_bg_cluster_tag       (0)
        , m_precomputed_marginal (_data.sizes(), 0)
        , m_data                 (&_data)
        , m_verbose              (verbose)
{
        fill(alpha.begin(), alpha.end(), 0.5);

        if (m_verbose >= 1) {
                flockfile(stderr);
                cerr << "Background gamma shape: " << parameters[0] << endl
                     << "Background gamma scale: " << parameters[1] << endl;
                fflush(stderr);
                funlockfile(stderr);
        }
}

default_background_t::default_background_t(const default_background_t& distribution)
        : component_model_t       (distribution)
        , alpha                  (distribution.alpha)
        , prior_distribution     (distribution.prior_distribution)
        , m_size                 (distribution.m_size)
        , m_bg_cluster_tag       (distribution.m_bg_cluster_tag)
        , m_precomputed_marginal (distribution.m_precomputed_marginal)
        , m_data                 (distribution.m_data)
        , m_verbose              (distribution.m_verbose)
{ }

default_background_t::~default_background_t() {
}

default_background_t*
default_background_t::clone() const {
        return new default_background_t(*this);
}

default_background_t&
default_background_t::operator=(const component_model_t& component_model)
{
        default_background_t tmp(
                static_cast<const default_background_t&>(component_model));
        swap(*this, tmp);
        return *this;
}

void
default_background_t::update(const string& msg_prefix)
{
        gradient_ascent();
        precompute_marginal();

        if (m_verbose >= 1) {
                flockfile(stderr);
                if (msg_prefix != "") {
                        cerr << msg_prefix << ": ";
                }
                cerr << "Background pseudocounts: ";
                for (size_t k = 0; k < m_size; k++) {
                        cerr << alpha[k] << " ";
                }
                cerr << endl;
                fflush(stderr);
                funlockfile(stderr);
        }
}

void
default_background_t::precompute_marginal()
{
        /* go through the data and precompute
         * lnbeta(n + alpha) - lnbeta(alpha) */
        for(size_t i = 0; i < data().size(); i++) {
                for(size_t j = 0; j < data()[i].size(); j++) {
                        m_precomputed_marginal[i][j] =
                                  mbeta_log(alpha, data()[i][j])
                                - mbeta_log(alpha);
                }
        }
}

double
default_background_t::gradient(
        const index_t& index, size_t k,
        double alpha_sum)
{
        double sum = accumulate(data()[index].begin(), data()[index].end(), 0.0);

        return boost::math::digamma(data()[index][k]+alpha[k]) - boost::math::digamma(sum+alpha_sum);
}

void
default_background_t::gradient(
        const index_t& index,
        double alpha_sum,
        vector<double>& result)
{
        for (size_t k = 0; k < m_size; k++) {
                result[k] += gradient(index, k, alpha_sum);
        }
}

void
default_background_t::gradient(vector<double>& result)
{
        double alpha_sum = accumulate(alpha.begin(), alpha.end(), 0.0);
        double n = 0.0;

        // likelihood
        for (size_t i = 0; i < data().size(); i++) {
                for (size_t j = 0; j < data()[i].size(); j++) {
                        index_t index(i, j);
                        if (cluster_assignments()[index] == m_bg_cluster_tag) {
                                gradient(index, alpha_sum, result);
                                n += 1.0;
                        }
                }
        }
        for (size_t k = 0; k < m_size; k++) {
                result[k] -= n*(boost::math::digamma(alpha[k]) - boost::math::digamma(alpha_sum));
        }
        // prior
        for (size_t k = 0; k < m_size; k++) {
                result[k] += boost::math::log_pdf_derivative(prior_distribution, alpha[k]);
        }
}

double
default_background_t::gradient_ascent(
        vector<double>& g,
        vector<double>& g_prev,
        vector<double>& epsilon,
        double eta,
        double min_alpha)
{
        double result = 0.0;

        /* save old gradient and compute the new one */
        g_prev = g; fill(g.begin(), g.end(), 0.0);
        gradient(g);

        for (size_t k = 0; k < m_size; k++) {
                double step = g[k] >= 0.0 ?
                        epsilon[k] : -epsilon[k];

                if (g[k] == 0.0)
                        continue;
                if (alpha[k] + step > 0.0) {
                        alpha[k] += step;
                        if (g_prev[k]*g[k] > 0.0) {
                                epsilon[k] *= 1.0+eta;
                        }
                        if (g_prev[k]*g[k] < 0.0) {
                                epsilon[k] *= 1.0-eta;
                        }
                }
                else {
                        alpha[k] = min_alpha;
                }
                result += abs(g[k]);
        }
        return result;
}

void
default_background_t::gradient_ascent()
{
        vector<double> g      (m_size, 1.0);
        vector<double> g_prev (m_size, 0.0);
        vector<double> epsilon(m_size, 1.0e-2);

        for (double sum = 1.0; sum > 0.1;) {
                sum = gradient_ascent(g, g_prev, epsilon);
        }
}

size_t
default_background_t::add(const range_t& range) {
        return range.length();
}

size_t
default_background_t::remove(const range_t& range) {
        return range.length();
}

size_t
default_background_t::count(const range_t& range) {
        return range.length();
}

/*
 *  p(y|x) = Beta(n(x) + n(y) + alpha) / Beta(n(x) + alpha)
 */
double default_background_t::predictive(const range_t& range) {
        return exp(log_predictive(range));
}

double default_background_t::predictive(const vector<range_t>& range_set) {
        return exp(log_predictive(range_set));
}

double default_background_t::log_predictive(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();
        double result = 0;

        for (size_t i = 0; i < length; i++) {
                const index_t index(sequence, position+i);

                /* counts contains the data count statistic
                 * and the pseudo counts alpha */
                result += m_precomputed_marginal[index];
        }

        return result;
}

double default_background_t::log_predictive(const vector<range_t>& range_set) {
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
                        result += m_precomputed_marginal[index];
                }
        }

        return result;
}

/*
 *  p(x) = Beta(n(x) + alpha) / Beta(alpha)
 */
double default_background_t::log_likelihood() const {
        double result = 0;

        /* counts contains the data count statistic
         * and the pseudo counts alpha */
        for(size_t i = 0; i < cluster_assignments().size(); i++) {
                for(size_t j = 0; j < cluster_assignments()[i].size(); j++) {
                        if (cluster_assignments()[i][j] == m_bg_cluster_tag) {
                                const index_t index(i, j);
                                result += m_precomputed_marginal[index];
                        }
                }
        }
        /* prior probability for the pseudocounts*/
        for (size_t k = 0; k < m_size; k++) {
                result += std::log(boost::math::pdf(prior_distribution, alpha[k]));
        }

        return result;
}

string
default_background_t::print_counts() const {
        return string();
}

void
default_background_t::set_bg_cluster_tag(cluster_tag_t bg_cluster_tag) {
        m_bg_cluster_tag = bg_cluster_tag;
}
