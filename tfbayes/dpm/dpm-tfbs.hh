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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/dpm/mixture-model.hh>
#include <tfbayes/dpm/data-tfbs.hh>
#include <tfbayes/dpm/dpm-tfbs-options.hh>
#include <tfbayes/dpm/dpm-tfbs-prior.hh>
#include <tfbayes/dpm/dpm-tfbs-state.hh>
#include <tfbayes/utility/linalg.h>

class dpm_tfbs_t : public mixture_model_t {
public:
         dpm_tfbs_t(const tfbs_options_t& options, const data_tfbs_t& data);
         dpm_tfbs_t(const dpm_tfbs_t& dpm);
        ~dpm_tfbs_t();

        virtual dpm_tfbs_t* clone() const;

        // auxiliary types
        ////////////////////////////////////////////////////////////////////////
        typedef mixture_state_t::const_iterator cm_iterator;
        typedef cluster_t::const_iterator cl_iterator;
        typedef data_tfbs_t::const_iterator da_iterator;

        // operators
        ////////////////////////////////////////////////////////////////////////
        friend std::ostream& operator<<(std::ostream& o, const dpm_tfbs_t& dpm);

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t mixture_components() const;
        size_t baseline_components() const;
        void   mixture_weights(const index_i& index, double log_weights[], cluster_tag_t tags[]);
        void   mixture_weights(const std::vector<range_t>& range_set, double log_weights[], cluster_tag_t cluster_tags[], bool include_background = true);
        void   update_map();
        void   update_graph(sequence_data_t<short> tfbs_start_positions);
        void   update_samples(size_t sampling_steps);
        double likelihood() const;
        double posterior() const;
        bool   valid_for_sampling(const index_i& index) const;

        const dpm_tfbs_state_t& state() const;
              dpm_tfbs_state_t& state();

              samples_t& samples();

        // test methods
        ////////////////////////////////////////////////////////////////////////
        void test();
        void test_background();
        void test_moves();
        void test_metropolis_hastings();

        // constants
        ////////////////////////////////////////////////////////////////////////
        static const size_t BG_LENGTH = 1;

protected:
        // baseline models
        std::vector<double> _baseline_weights;
        std::vector<baseline_tag_t> _baseline_tags;

        // data and clusters
        const data_tfbs_t& _data;
        dpm_tfbs_state_t _state;

        // parameters
        const double _lambda;
        const double _lambda_log;
        const double _lambda_inv_log;
        const size_t _tfbs_length;

        // samples
        samples_t _samples;

        // process priors
        dpm_tfbs_prior_t* _process_prior;

        // whether or not to construct a graph from the samples
        bool _construct_graph;
};

#endif /* DPM_TFBS_HH */
