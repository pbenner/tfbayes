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

#ifndef DPM_TFBS_HH
#define DPM_TFBS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <clustermanager.hh>
#include <component-model.hh>
#include <data-tfbs.hh>
#include <mixture-model.hh>

class DPM_TFBS : public DPM {
public:
         DPM_TFBS(double alpha, double lambda, size_t tfbs_length, const DataTFBS& data);
        ~DPM_TFBS();

        DPM_TFBS* clone() const;

        // operators
        ////////////////////////////////////////////////////////////////////////
              Cluster& operator[](cluster_tag_t c)       { return _cluster_manager[c]; }
        const Cluster& operator[](cluster_tag_t c) const { return _cluster_manager[c]; }

        friend std::ostream& operator<<(std::ostream& o, const DPM_TFBS& dpm);

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
        const DataTFBS& data() const;
        const ClusterManager& cluster_manager() const;

        // constants
        ////////////////////////////////////////////////////////////////////////
        static const size_t BG_LENGTH   = 1;
               const size_t TFBS_LENGTH;
        static const size_t NUCLEOTIDES = 4;

private:
        // priors
        gsl_matrix* bg_alpha;
        gsl_matrix* tfbs_alpha;

        // data and clusters
        const DataTFBS& _data;
        sequence_data_t<cluster_tag_t> _cluster_assignments;
        ClusterManager   _cluster_manager;

        // tags of special clusters
        cluster_tag_t bg_cluster_tag;

        // parameters
        double alpha;
        double lambda;

        // record start positions of tfbs
        sequence_data_t<short> _tfbs_start_positions;

        // posterior distribution
        posterior_t _posterior;

        // keep track of the number of transcription factor binding sites
        size_t num_tfbs;

        // standard priors
        static gsl_matrix* init_alpha(size_t length) {
                gsl_matrix* bg_alpha = gsl_matrix_alloc(length, DPM_TFBS::NUCLEOTIDES);

                // initialize prior for the background model
                for (size_t i = 0; i < length; i++) {
                        for (size_t j = 0; j < DPM_TFBS::NUCLEOTIDES; j++) {
                                gsl_matrix_set(bg_alpha, i, j, 1);
                        }
                }
                return bg_alpha;
        }
};

#endif /* DPM_TFBS_HH */
