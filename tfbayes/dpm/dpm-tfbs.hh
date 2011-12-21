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
#include <component-model.hh>
#include <data-tfbs.hh>
#include <mixture-model.hh>
#include <dpm-tfbs-state.hh>

#include <tfbayes/linalg.h>

typedef struct {
        size_t tfbs_length;
        double alpha;
        double discount;
        double lambda;
        size_t context;
        std::string process_prior;
        std::string background_model;
        std::vector<double> baseline_weights;
        std::vector<std::matrix<double> > baseline_priors;
} tfbs_options_t;

class DpmTfbs : public DPM {
public:
         DpmTfbs(const tfbs_options_t& options, const data_tfbs_t& data);
        ~DpmTfbs();

        DpmTfbs* clone() const;

        // auxiliary types
        ////////////////////////////////////////////////////////////////////////
        typedef mixture_state_t::const_iterator cm_iterator;
        typedef Cluster::const_iterator cl_iterator;
        typedef data_tfbs_t::const_iterator da_iterator;

        // operators
        ////////////////////////////////////////////////////////////////////////
              Cluster& operator[](cluster_tag_t c)       { return _state[c]; }
        const Cluster& operator[](cluster_tag_t c) const { return _state[c]; }

        friend std::ostream& operator<<(std::ostream& o, const DpmTfbs& dpm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t mixture_components() const;
        size_t baseline_components() const;
        void   mixture_weights(const index_i& index, double log_weights[], cluster_tag_t tags[]);
        void   add(const index_i& index, cluster_tag_t tag);
        void   remove(const index_i& index, cluster_tag_t tag);
        void   update_graph(sequence_data_t<short> tfbs_start_positions);
        void   update_hypergraph(sequence_data_t<short> tfbs_start_positions);
        void   update_posterior(size_t sampling_steps);
        double likelihood() const;
        bool   valid_for_sampling(const index_i& index) const;
        posterior_t& posterior();
        const data_tfbs_t& data() const;
        const mixture_state_t& state() const;
        const Graph& graph() const;

        // test methods
        ////////////////////////////////////////////////////////////////////////
        void test();
        void test_background();
        void test_moves();
        void test_metropolis_hastings();

        // metropolis hastings
        void metropolis_hastings();
        bool proposal(Cluster& cluster);
        bool move_left(Cluster& cluster);
        bool move_right(Cluster& cluster);

        // constants
        ////////////////////////////////////////////////////////////////////////
        static const size_t BG_LENGTH   = 1;
               const size_t TFBS_LENGTH;
        static const size_t ALPHABET_SIZE = 4;

private:
        // baseline models
        std::vector<double> _baseline_weights;
        std::vector<model_tag_t> _model_tags;

        // data and clusters
        const data_tfbs_t& _data;
        dpm_tfbs_state_t _state;

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

        // process priors
        double py_prior(Cluster& cluster);
        double uniform_prior(Cluster& cluster);
        double poppe_prior(Cluster& cluster);
        typedef double (DpmTfbs::*prior_fn)(Cluster& cluster);
        prior_fn _process_prior;

        // standard priors
        static std::matrix<double> init_alpha(size_t length);
};

#endif /* DPM_TFBS_HH */
