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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_matrix.h>

#include <data-gaussian.hh>
#include <dpm.hh>
#include <distribution.hh>

class DPM_Gaussian : public DPM {
public:
         DPM_Gaussian(gsl_matrix* _cov,
                      gsl_matrix* _cov_0,
                      gsl_vector* _mu_0,
                      const Data_Gaussian& data);
        ~DPM_Gaussian();

        DPM_Gaussian* clone() const;

        // type definitions
        ////////////////////////////////////////////////////////////////////////
        typedef std::vector<std::vector<bool>   > tfbs_start_positions_t;

        // operators
        ////////////////////////////////////////////////////////////////////////
              Cluster& operator[](cluster_tag_t c)       { return _cluster_manager[c]; }
        const Cluster& operator[](cluster_tag_t c) const { return _cluster_manager[c]; }

        friend std::ostream& operator<<(std::ostream& o, const DPM_Gaussian& dpm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t mixture_components() const;
        void   mixture_weights(const index_t& index, double weights[], cluster_tag_t tags[]);
        void   add(const index_t& index, cluster_tag_t tag);
        void   remove(const index_t& index, cluster_tag_t tag);
        void   update_posterior(size_t sampling_steps);
        double likelihood() const;
        bool   valid_for_sampling(const index_t& index) const;
        const posterior_t& posterior() const;
        const ClusterManager& cluster_manager() const;

private:
        // likelihood parameters
        gsl_matrix* cov;
        gsl_matrix* cov_inv;

        // prior parameters
        gsl_vector* mu_0;
        gsl_matrix* cov_0;
        gsl_matrix* cov_inv_0;

        // data and clusters
        const Data_Gaussian& _data;
        data_t<cluster_tag_t> _cluster_assignments;
        ClusterManager _cluster_manager;

        // parameters
        double alpha;

        // posterior distribution
        posterior_t _posterior;
};

#endif /* DPM_GAUSSIAN_HH */
