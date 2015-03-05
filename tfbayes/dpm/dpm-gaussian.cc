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

#include <limits>

#include <tfbayes/dpm/dpm-gaussian.hh>
#include <tfbayes/utility/logarithmetic.hh>

using namespace std;

dpm_gaussian_t::dpm_gaussian_t(
        double alpha,
        const matrix<double>& Sigma,
        const matrix<double>& Sigma_0,
        const vector<double>& mu_0,
        const data_gaussian_t& data)
        : gibbs_state_t(data_t<cluster_tag_t>(data.elements(), -1)),
          _data(&data),
          // strength parameter for the dirichlet process
          _alpha(alpha)
{
        _baseline_tag = gibbs_state_t::add_baseline_model(new bivariate_normal_t(Sigma, Sigma_0, mu_0, data));
        cluster_tag_t cluster_tag = gibbs_state_t::get_free_cluster(_baseline_tag).cluster_tag();
        for (indexer_t::const_iterator it = _data->begin();
             it != _data->end(); it++) {
                gibbs_state_t::operator[](cluster_tag).add_observations(range_t(**it,1));
        }
}

dpm_gaussian_t::dpm_gaussian_t(const dpm_gaussian_t& dpm)
        : gibbs_state_t(dpm),
          _baseline_tag(dpm._baseline_tag),
          _data        (dpm._data),
          _alpha       (dpm._alpha)
{
}

dpm_gaussian_t::~dpm_gaussian_t() {
}

void swap(dpm_gaussian_t& first, dpm_gaussian_t& second)
{
        swap(static_cast<gibbs_state_t&>(first),
             static_cast<gibbs_state_t&>(second));
        swap(first._baseline_tag, second._baseline_tag);
        swap(first._data,         second._data);
        swap(first._alpha,        second._alpha);
}

dpm_gaussian_t*
dpm_gaussian_t::clone() const {
        return new dpm_gaussian_t(*this);
}

dpm_gaussian_t&
dpm_gaussian_t::operator=(const mixture_model_t& mixture_model)
{
        dpm_gaussian_t tmp(static_cast<const dpm_gaussian_t&>(mixture_model));
        swap(*this, tmp);
        return *this;
}

void
dpm_gaussian_t::add(const range_t& range, cluster_tag_t cluster_tag)
{
        gibbs_state_t::operator[](cluster_tag).add_observations(range);
}

void
dpm_gaussian_t::remove(const range_t& range)
{
        cluster_tag_t cluster_tag = state()[range.index()];

        gibbs_state_t::operator[](cluster_tag).remove_observations(range);
}

size_t
dpm_gaussian_t::mixture_components() const
{
        return gibbs_state_t::size();
}

void
dpm_gaussian_t::mixture_weights(const range_t& range, double log_weights[], cluster_tag_t cluster_tags[])
{
        size_t components = mixture_components();
        double sum        = -numeric_limits<double>::infinity();
        double N          = _data->elements() - 1;

        cluster_tag_t i = 0;
        for (mixture_state_t::const_iterator it = gibbs_state_t::begin(); it != gibbs_state_t::end(); it++) {
                cluster_t& cluster = **it;
                cluster_tags[i] = cluster.cluster_tag();
                double num_elements = (double)cluster.size();
                // normalization constant
                sum = logadd(sum, log(num_elements/(_alpha + N)) + cluster.model().log_predictive(range));
                log_weights[i] = sum;
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        cluster_tags[components] = gibbs_state_t::get_free_cluster(_baseline_tag).cluster_tag();
        sum = logadd(sum, log(_alpha/(_alpha + N)) + gibbs_state_t::operator[](cluster_tags[components]).model().log_predictive(range));
        log_weights[components] = sum;
}

matrix<double>
dpm_gaussian_t::means() const {
        matrix<double> means;

        for (mixture_state_t::const_iterator it = gibbs_state_t::begin();
             it != gibbs_state_t::end(); it++) {
                cluster_t& cluster = **it;
                bivariate_normal_t& bg = static_cast<bivariate_normal_t&>(cluster.model());
                means.push_back(bg.mean());
        }

        return means;
}

double
dpm_gaussian_t::likelihood() const {
        double result = 0;

        for (mixture_state_t::const_iterator it = gibbs_state_t::begin();
             it != gibbs_state_t::end(); it++) {
                cluster_t& cluster = **it;
                result += cluster.model().log_likelihood();
        }
        return result;
}

double
dpm_gaussian_t::posterior() const {
        return 0;
}

const mixture_state_t&
dpm_gaussian_t::state() const {
        return *this;
}

mixture_state_t&
dpm_gaussian_t::state() {
        return *this;
}

// misc methods
////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const dpm_gaussian_t& dpm)
{
        return o;
}
