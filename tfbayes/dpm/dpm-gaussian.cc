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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <tfbayes/dpm/dpm-gaussian.hh>
#include <tfbayes/utility/logarithmetic.h>

using namespace std;

dpm_gaussian_t::dpm_gaussian_t(
        double alpha,
        gsl_matrix* Sigma,
        gsl_matrix* Sigma_0,
        gsl_vector* mu_0,
        const data_gaussian_t& data)
        : gibbs_state_t(data_t<cluster_tag_t>(data.elements(), -1)),
          _data(data),
          // strength parameter for the dirichlet process
          alpha(alpha)
{
        _baseline_tag = "baseline";
        gibbs_state_t::add_baseline_model(new bivariate_normal_t(Sigma, Sigma_0, mu_0, data), _baseline_tag);
        cluster_tag_t cluster_tag = gibbs_state_t::get_free_cluster(_baseline_tag).cluster_tag();
        for (data_gaussian_t::const_iterator it = _data.begin();
             it != _data.end(); it++) {
                gibbs_state_t::operator[](cluster_tag).add_observations(range_t(**it,1));
        }
}

dpm_gaussian_t::~dpm_gaussian_t() {
}

dpm_gaussian_t*
dpm_gaussian_t::clone() const {
        return new dpm_gaussian_t(*this);
}

bool
dpm_gaussian_t::valid_for_sampling(const index_i& index) const
{
        return true;
}

void
dpm_gaussian_t::add(const index_i& index, cluster_tag_t cluster_tag)
{
        gibbs_state_t::operator[](cluster_tag).add_observations(range_t(index,1));
}

void
dpm_gaussian_t::remove(const index_i& index, cluster_tag_t cluster_tag)
{
        gibbs_state_t::operator[](cluster_tag).remove_observations(range_t(index,1));
}

size_t
dpm_gaussian_t::mixture_components() const
{
        return gibbs_state_t::size();
}

void
dpm_gaussian_t::mixture_weights(const index_i& index, double log_weights[], cluster_tag_t cluster_tags[])
{
        size_t components = mixture_components();
        double sum        = -HUGE_VAL;
        double N          = _data.elements() - 1;
        range_t range(index, 1);

        cluster_tag_t i = 0;
        for (mixture_state_t::const_iterator it = gibbs_state_t::begin(); it != gibbs_state_t::end(); it++) {
                cluster_t& cluster = **it;
                cluster_tags[i] = cluster.cluster_tag();
                double num_elements = (double)cluster.size();
                // normalization constant
                sum = logadd(sum, log(num_elements/(alpha + N)) + cluster.model().log_predictive(range));
                log_weights[i] = sum;
                i++;
        }
        ////////////////////////////////////////////////////////////////////////
        // add the tag of a new class and compute their weight
        cluster_tags[components] = gibbs_state_t::get_free_cluster(_baseline_tag).cluster_tag();
        sum = logadd(sum, log(alpha/(alpha + N)) + gibbs_state_t::operator[](cluster_tags[components]).model().log_predictive(range));
        log_weights[components] = sum;
}

gsl_matrix*
dpm_gaussian_t::means() const {
        if (mixture_components() == 0) {
                return NULL;
        }

        gsl_matrix* means = gsl_matrix_alloc(mixture_components(), 2);

        size_t i = 0;
        for (mixture_state_t::const_iterator it = gibbs_state_t::begin();
             it != gibbs_state_t::end(); it++) {
                cluster_t& cluster = **it;
                bivariate_normal_t& bg = static_cast<bivariate_normal_t&>(cluster.model());
                gsl_matrix_set(means, i, 0, gsl_vector_get(bg.mean(), 0));
                gsl_matrix_set(means, i, 1, gsl_vector_get(bg.mean(), 1));
                i++;
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

void
dpm_gaussian_t::update_samples(size_t sampling_steps) {
}

samples_t&
dpm_gaussian_t::samples() {
        return _samples;
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
