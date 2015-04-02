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

#include <algorithm>    // std::accumulate
#include <cmath>        // std::exp, std::log
#include <string>
#include <sstream>
#include <vector>

#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/utility/multinomial-beta.hh>
#include <tfbayes/utility/normalize.hh>
#include <tfbayes/utility/random.hh>

#include <boost/math/distributions/gamma_extra.hpp>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

default_background_t::default_background_t(
        const vector<double>& parameters,
        const vector<double>& weights,
        const sequence_data_t<data_tfbs_t::code_t>& data,
        const sequence_data_t<cluster_tag_t>& cluster_assignments,
        thread_pool_t& thread_pool,
        const string& cachefile,
        boost::optional<const alignment_set_t<>&> alignment_set,
        size_t verbose)
        : component_model_t       ({"background", 1}, cluster_assignments)
        , m_alpha                 (weights.size(), counts_t())
        , m_weights               (normalize(weights))
        , m_prior_distribution    (parameters[0], parameters[1])
        , m_size1                 (weights.size())
        , m_size2                 (data_tfbs_t::alphabet_size)
        , m_bg_cluster_tag        (0)
        , m_log_likelihood        (0.0)
        , m_marginal_probability  (data.sizes(),  0)
        , m_component_assignments (data.sizes(), -1)
        , m_data                  (&data)
        , m_verbose               (verbose)
{
        boost::random::mt19937 gen; seed_rng(gen);
        // distribution for generating noise
        boost::random::normal_distribution<> dist(0.0, 0.001);
        // if no weights are given, assume the model should consist of
        // a single component
        if (m_size1 == 0) {
                m_size1 = 1;
                m_alpha  .push_back(counts_t());
                m_weights.push_back(1.0);
        }
        // initialize parameters differently for each component
        for (size_t i = 0; i < m_size1; i++) {
                for (size_t j = 0; j < m_size2; j++) {
                        m_alpha[i][j] = 0.5 + dist(gen);
                }
        }

        if (m_verbose >= 1) {
                flockfile(stderr);
                cerr << "Background gamma shape: " << parameters[0] << endl
                     << "Background gamma scale: " << parameters[1] << endl
                     << "Background components : " << m_size1       << endl;
                fflush(stderr);
                funlockfile(stderr);
        }
        // initialize gradients
        m_g       = matrix<double>(m_size1, m_size2, 1.0);
        m_g_prev  = matrix<double>(m_size1, m_size2, 0.0);
        m_epsilon = matrix<double>(m_size1, m_size2, 1.0e-2);
        m_n       = vector<double>(m_size1, 0.0);
}

default_background_t::default_background_t(const default_background_t& distribution)
        : component_model_t       (distribution)
        , m_alpha                 (distribution.m_alpha)
        , m_weights               (distribution.m_weights)
        , m_prior_distribution    (distribution.m_prior_distribution)
        , m_size1                 (distribution.m_size1)
        , m_size2                 (distribution.m_size2)
        , m_bg_cluster_tag        (distribution.m_bg_cluster_tag)
        , m_log_likelihood        (distribution.m_log_likelihood)
        , m_marginal_probability  (distribution.m_marginal_probability)
        , m_component_assignments (distribution.m_component_assignments)
        , m_data                  (distribution.m_data)
        , m_verbose               (distribution.m_verbose)
        , m_g                     (distribution.m_g)
        , m_g_prev                (distribution.m_g_prev)
        , m_epsilon               (distribution.m_epsilon)
        , m_n                     (distribution.m_n)
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
default_background_t::compute_marginal()
{
        /* also update the log likelihood */
        m_log_likelihood = 0.0;

        /* go through the data and precompute
         * lnbeta(n + alpha) - lnbeta(alpha) */
        for(size_t i = 0; i < data().size(); i++) {
                for(size_t j = 0; j < data()[i].size(); j++) {
                        const index_t index(i,j);
                        /* get mixture component */
                        const size_t k = m_component_assignments[index];
                        /* recompute marginal at this position */
                        m_marginal_probability[index] =
                                  mbeta_log(m_alpha[k], data()[index])
                                - mbeta_log(m_alpha[k]);
                        /* if this position is assigned to the
                         * background, update the log likelihood */
                        if (cluster_assignments()[index] == m_bg_cluster_tag) {
                                m_log_likelihood += m_marginal_probability[index];
                        }
                }
        }
}

void
default_background_t::gradient(
        const index_t& index,
        size_t i, size_t j,
        double alpha_sum)
{
        const double sum = accumulate(data()[index].begin(), data()[index].end(), 0.0);

        m_g[i][j] += boost::math::digamma(data()[index][j]+m_alpha[i][j])
                - boost::math::digamma(sum+alpha_sum);
}

void
default_background_t::gradient(
        const index_t& index,
        size_t i,
        double alpha_sum)
{
        for (size_t j = 0; j < m_size2; j++) {
                gradient(index, i, j, alpha_sum);
        }
}

void
default_background_t::gradient()
{
        // precompute pseudocount sums for each component
        vector<double> alpha_sum;
        for (size_t i = 0; i < m_size1; i++) {
                alpha_sum.push_back(accumulate(m_alpha[i].begin(), m_alpha[i].end(), 0.0));
        };

        // likelihood
        for (size_t i = 0; i < data().size(); i++) {
                for (size_t j = 0; j < data()[i].size(); j++) {
                        index_t index(i, j);
                        /* the gradient should consider only those
                         * positions that are currently assigned to the
                         * background */
                        if (cluster_assignments()[index] == m_bg_cluster_tag) {
                                const ssize_t k = m_component_assignments[index];
                                assert(k != -1);

                                gradient(index, k, alpha_sum[k]);
                        }
                }
        }
        for (size_t i = 0; i < m_size1; i++) {
                for (size_t j = 0; j < m_size2; j++) {
                        m_g[i][j] -= m_n[i]*(boost::math::digamma(m_alpha[i][j]) - boost::math::digamma(alpha_sum[i]));
                }
        }
        // prior
        for (size_t i = 0; i < m_size1; i++) {
                for (size_t j = 0; j < m_size2; j++) {
                        m_g[i][j] += boost::math::log_pdf_derivative(m_prior_distribution, m_alpha[i][j]);
                }
        }
}

double
default_background_t::gradient_ascent_loop(
        double eta,
        double min_alpha)
{
        double result = 0.0;

        /* save old gradient and compute the new one */
        m_g_prev = m_g;
        /* reset gradient */
        for (size_t i = 0; i < m_size1; i++) {
                fill(m_g[i].begin(), m_g[i].end(), 0.0);
        }
        /* recompute gradient */
        gradient();

        /* recompute step size and
         * update alpha pseudocounts*/
        for (size_t i = 0; i < m_size1; i++) {
                for (size_t j = 0; j < m_size2; j++) {
                        double step = m_g[i][j] >= 0.0 ?
                                m_epsilon[i][j] : -m_epsilon[i][j];

                        if (m_g[i][j] == 0.0)
                                continue;
                        if (m_alpha[i][j] + step > 0.0) {
                                m_alpha[i][j] += step;
                                if (m_g_prev[i][j]*m_g[i][j] > 0.0) {
                                        m_epsilon[i][j] *= 1.0+eta;
                                }
                                if (m_g_prev[i][j]*m_g[i][j] < 0.0) {
                                        m_epsilon[i][j] *= 1.0-eta;
                                }
                        }
                        else {
                                m_alpha[i][j] = min_alpha;
                        }
                        result += abs(m_g[i][j]);
                }
        }
        return result;
}

bool
default_background_t::gradient_ascent()
{
        bool optimized = false;

        for (double sum = 1.0; sum > 0.1;) {
                sum = gradient_ascent_loop();
                /* record of parameters changed */
                if (sum > 0.1) {
                        optimized = true;
                }
        }
        return optimized;
}

ssize_t
default_background_t::max_component(const index_t& index) const
{
        vector<double> result(m_size1, 0.0);

        /* find component with highest predictive value */
        for (size_t i = 0; i < m_size1; i++) {
                /* counts contains the data count statistic
                 * and the pseudo counts alpha */
                result[i] = mbeta_log(m_alpha[i], data()[index])
                          - mbeta_log(m_alpha[i])
                          + log(m_weights[i]);
        }
        return distance(result.begin(), max_element(result.begin(), result.end()));
}

bool
default_background_t::compute_component_assignments_loop()
{
        bool optimized = false;
        /* optimize assignments */
        for (size_t i = 0; i < data().size(); i++) {
                for (size_t j = 0; j < data()[i].size(); j++) {
                        index_t index(i, j);
                        /* update count statistics */
                        if (cluster_assignments()[index] == m_bg_cluster_tag && 
                            m_component_assignments[index] != -1) {
                                m_n[m_component_assignments[index]] -= 1.0;
                        }
                        /* get best assignment */
                        ssize_t k = max_component(index);
                        /* check if assigment changed */
                        if (k != m_component_assignments[index]) {
                                optimized = true;
                        }
                        /* save assignemnt */
                        m_component_assignments[index] = k;
                        /* update count statistics */
                        if (cluster_assignments()[index] == m_bg_cluster_tag) {
                                m_n[k] += 1.0;
                        }
                }
        }
        return optimized;
}

bool
default_background_t::compute_component_assignments()
{
        bool optimized = true;
        bool result    = false;
        /* iterate until reaching a fixed point */
        for (size_t i = 0; optimized; i++) {
                optimized = compute_component_assignments_loop();
                result |= optimized;
        }
        return result;
}

void
default_background_t::update(const string& msg_prefix)
{
        bool optimized;

        do {
                optimized = false;
                /* given the pseucodounts, recompute component
                 * assignments */
                optimized |= compute_component_assignments();
                /* optimized pseudocounts */
                optimized |= gradient_ascent();
        } while (optimized);

        /* recompute marginal probabilities and the log likelihood */
        compute_marginal();

        if (m_verbose >= 1) {
                flockfile(stderr);
                if (msg_prefix != "") {
                        cerr << msg_prefix << ": ";
                }
                cerr << "Background pseudocounts: "
                     << endl;
                cerr << print_pseudocounts()
                     << endl;
                fflush(stderr);
                funlockfile(stderr);
        }
}

size_t
default_background_t::add(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();

        for (size_t i = 0; i < length; i++) {
                const index_t index(sequence, position+i);
                m_log_likelihood += m_marginal_probability[index];
                /* update count statistics */
                if (m_component_assignments[index] != -1) {
                        m_n[m_component_assignments[index]] += 1.0;
                }
        }

        return length;
}

size_t
default_background_t::remove(const range_t& range) {
        const size_t sequence = range.index()[0];
        const size_t position = range.index()[1];
        const size_t length   = range.length();

        for (size_t i = 0; i < length; i++) {
                const index_t index(sequence, position+i);
                m_log_likelihood -= m_marginal_probability[index];
                /* update count statistics */
                if (m_component_assignments[index] != -1) {
                        m_n[m_component_assignments[index]] -= 1.0;
                }
        }

        return length;
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
                result += m_marginal_probability[index];
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
                        result += m_marginal_probability[index];
                }
        }

        return result;
}

/*
 *  p(x) = Beta(n(x) + alpha) / Beta(alpha)
 */
double default_background_t::log_likelihood() const {
        double result = m_log_likelihood;

        /* pseudocount probability */
        for (size_t i = 0; i < m_size1; i++) {
                for (size_t j = 0; j < m_size2; j++) {
                        result += std::log(boost::math::pdf(m_prior_distribution, m_alpha[i][j]));
                }
        }

        return result;
}

string
default_background_t::print_pseudocounts() const {
        stringstream ss;

        for (size_t i = 0; i < m_size1; i++) {
                if (i != 0) ss << endl;
                ss << " -> ";
                for (size_t j = 0; j < m_size2; j++) {
                        ss << boost::format("%7.4f ") % m_alpha[i][j];
                }
                ss << boost::format("(%0.2f%%)") % (m_n[i]/accumulate(m_n.begin(), m_n.end(), 0.0));
        }
        return ss.str();
}

string
default_background_t::print_counts() const {
        return string();
}

void
default_background_t::set_bg_cluster_tag(cluster_tag_t bg_cluster_tag) {
        m_bg_cluster_tag = bg_cluster_tag;
}
