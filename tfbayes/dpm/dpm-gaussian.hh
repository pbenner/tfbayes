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

#ifndef __TFBAYES_DPM_DPM_GAUSSIAN_HH__
#define __TFBAYES_DPM_DPM_GAUSSIAN_HH__

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
                        const std::matrix<double>& Sigma,
                        const std::matrix<double>& Sigma_0,
                        const std::vector<double>& mu_0,
                        const data_gaussian_t& data);
         dpm_gaussian_t(const dpm_gaussian_t& dpm);
        ~dpm_gaussian_t();

        virtual dpm_gaussian_t* clone() const;

        friend void swap(dpm_gaussian_t& first, dpm_gaussian_t& second);

        // operators
        ////////////////////////////////////////////////////////////////////////
        friend std::ostream& operator<<(std::ostream& o, const dpm_gaussian_t& dpm);

        virtual dpm_gaussian_t& operator=(const mixture_model_t& mixture_model);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t baseline_components() const { return 1; }
        size_t mixture_components() const;
        void   mixture_weights(const range_t& range, double log_weights[], cluster_tag_t tags[]);
        void   add   (const range_t& range, cluster_tag_t tag);
        void   remove(const range_t& range, cluster_tag_t tag);
        void   remove(const index_i& index, cluster_tag_t tag);
        double likelihood() const;
        double posterior() const;
        std::matrix<double> means() const;

        const mixture_state_t& state() const;
              mixture_state_t& state();

private:
        baseline_tag_t _baseline_tag;

        // data and clusters
        const data_gaussian_t* _data;

        // parameters
        double _alpha;
};

#endif /* __TFBAYES_DPM_DPM_GAUSSIAN_HH__ */
