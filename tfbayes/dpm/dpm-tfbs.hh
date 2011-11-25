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

#include <graph.hh>
#include <clustermanager.hh>
#include <component-model.hh>
#include <data-tfbs.hh>
#include <mixture-model.hh>

class DPM_TFBS : public DPM {
public:
         DPM_TFBS(double alpha, double discount, double lambda, size_t tfbs_length,
                  const data_tfbs_t& data, const data_tfbs_t& data_comp,
                  std::vector<double> baseline_weights, gsl_matrix *baseline_priors[],
                  std::string process_prior_name = "pitman-yor process");
        ~DPM_TFBS();

        DPM_TFBS* clone() const;

        // operators
        ////////////////////////////////////////////////////////////////////////
              Cluster& operator[](cluster_tag_t c)       { return _clustermanager[c]; }
        const Cluster& operator[](cluster_tag_t c) const { return _clustermanager[c]; }

        friend std::ostream& operator<<(std::ostream& o, const DPM_TFBS& dpm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t mixture_components() const;
        size_t baseline_components() const;
        void   mixture_weights(const index_t& index, double log_weights[], cluster_tag_t tags[]);
        void   add(const index_t& index, cluster_tag_t tag);
        void   remove(const index_t& index, cluster_tag_t tag);
        void   update_graph(sequence_data_t<short> tfbs_start_positions);
        void   update_hypergraph(sequence_data_t<short> tfbs_start_positions);
        void   update_posterior(size_t sampling_steps);
        double likelihood() const;
        bool   valid_for_sampling(const index_t& index) const;
        posterior_t& posterior();
        const data_tfbs_t& data() const;
        const ClusterManager& clustermanager() const;
        const Graph& graph() const;
        void test();

        // constants
        ////////////////////////////////////////////////////////////////////////
        static const size_t BG_LENGTH   = 1;
               const size_t TFBS_LENGTH;
        static const size_t PARSMM_DEPTH = 3;
        static const size_t ALPHABET_SIZE = 4;

private:
        // baseline models
        std::vector<double> _baseline_weights;
        std::vector<model_tag_t> _model_tags;

        // data and clusters
        const data_tfbs_t& _data;
        const data_tfbs_t& _data_comp;
        sequence_data_t<cluster_tag_t> _cluster_assignments;
        ClusterManager _clustermanager;

        // tags of special clusters
        cluster_tag_t bg_cluster_tag;

        // parameters
        const double alpha;
        const double alpha_log;
        const double discount;
        const double discount_log;
        const double lambda;
        const double lambda_log;
        const double lambda_inv_log;

        // record start positions of tfbs
        sequence_data_t<short> _tfbs_start_positions;

        // posterior distribution
        posterior_t _posterior;
        Graph _tfbs_graph;

        // keep track of the number of transcription factor binding sites
        size_t num_tfbs;

        // process priors
        double py_prior(Cluster& cluster);
        double uniform_prior(Cluster& cluster);
        double poppe_prior(Cluster& cluster);
        typedef double (DPM_TFBS::*prior_fn)(Cluster& cluster);
        prior_fn _process_prior;

        // standard priors
        static gsl_matrix* init_alpha(size_t length) {
                gsl_matrix* alpha = gsl_matrix_alloc(length, ALPHABET_SIZE);

                // initialize prior for the background model
                for (size_t i = 0; i < length; i++) {
                        for (size_t j = 0; j < ALPHABET_SIZE; j++) {
                                gsl_matrix_set(alpha, i, j, 1);
                        }
                }
                return alpha;
        }
};

#endif /* DPM_TFBS_HH */
