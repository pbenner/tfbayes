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

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <clustermanager.hh>
#include <dpm.hh>
#include <distribution.hh>

class DPM_GAUSSIAN : public DPM {
public:
         DPM_GAUSSIAN(gsl_matrix* _cov,
                      gsl_matrix* _cov_0,
                      gsl_vector* _mu_0,
                      const Data& data);
        ~DPM_GAUSSIAN();

        DPM_GAUSSIAN* clone() const;

        // type definitions
        ////////////////////////////////////////////////////////////////////////
        typedef std::vector<std::vector<bool>   > tfbs_start_positions_t;

        // operators
        ////////////////////////////////////////////////////////////////////////
              Cluster& operator[](cluster_tag_t c)       { return _cluster_manager[c]; }
        const Cluster& operator[](cluster_tag_t c) const { return _cluster_manager[c]; }

        friend std::ostream& operator<<(std::ostream& o, const DPM_GAUSSIAN& dpm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t mixture_components() const;
        void   mixture_weights(const word_t& word, double weights[], cluster_tag_t tags[]);
        void   add_word(const word_t& word, cluster_tag_t tag);
        void   remove_word(const word_t& word, cluster_tag_t tag);
        void   update_posterior(size_t sampling_steps);
        double likelihood() const;
        bool   valid_for_sampling(const element_t& element, const word_t& word) const;
        const posterior_t& posterior() const;
        const Data& data() const;
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
        const Data& _data;
        ClusterManager _cluster_manager;

        // parameters
        double alpha;

        // posterior distribution
        posterior_t _posterior;
};

#endif /* DPM_GAUSSIAN_HH */
