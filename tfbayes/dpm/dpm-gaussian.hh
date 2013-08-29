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

#ifndef DPM_GAUSSIAN_HH
#define DPM_GAUSSIAN_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_matrix.h>

#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/dpm/data-gaussian.hh>
#include <tfbayes/dpm/mixture-model.hh>
#include <tfbayes/dpm/state.hh>

class dpm_gaussian_t : public mixture_model_t, public gibbs_state_t {
public:
         dpm_gaussian_t(double alpha,
                      gsl_matrix* _Sigma,
                      gsl_matrix* _Sigma_0,
                      gsl_vector* _mu_0,
                      const data_gaussian_t& data);
        ~dpm_gaussian_t();

        dpm_gaussian_t* clone() const;

        // operators
        ////////////////////////////////////////////////////////////////////////
        friend std::ostream& operator<<(std::ostream& o, const dpm_gaussian_t& dpm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t baseline_components() const { return 1; }
        size_t mixture_components() const;
        void   mixture_weights(const index_i& index, double log_weights[], cluster_tag_t tags[]);
        void   add(const index_i& index, cluster_tag_t tag);
        void   remove(const index_i& index, cluster_tag_t tag);
        void   update_samples(size_t sampling_steps);
        double likelihood() const;
        double posterior() const;
        bool   valid_for_sampling(const index_i& index) const;
        gsl_matrix* means() const;

        const mixture_state_t& state() const;
              mixture_state_t& state();
              samples_t& samples();

        void print(std::ostream& o) const {
                o << "(" << gibbs_state_t::size() << "): ";
                for (mixture_state_t::const_iterator it = gibbs_state_t::begin();
                     it != gibbs_state_t::end(); it++) {
                        o << **it << " ";
                }
        }

private:
        baseline_tag_t _baseline_tag;

        // likelihood parameters
        gsl_matrix* cov;
        gsl_matrix* cov_inv;

        // prior parameters
        gsl_vector* mu_0;
        gsl_matrix* cov_0;
        gsl_matrix* cov_inv_0;

        // data and clusters
        const data_gaussian_t& _data;

        // parameters
        const double alpha;

        // samples
        samples_t _samples;
};

#endif /* DPM_GAUSSIAN_HH */
