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

#ifndef __TFBAYES_DPM_DPM_TFBS_HH__
#define __TFBAYES_DPM_DPM_TFBS_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <limits>

#include <boost/optional.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/dpm/mixture-model.hh>
#include <tfbayes/dpm/data-tfbs.hh>
#include <tfbayes/dpm/dpm-tfbs-options.hh>
#include <tfbayes/dpm/dpm-tfbs-prior.hh>
#include <tfbayes/dpm/dpm-tfbs-state.hh>
#include <tfbayes/dpm/dpm-sampling-history.hh>
#include <tfbayes/utility/linalg.hh>

class dpm_tfbs_t : public mixture_model_t {
public:
         dpm_tfbs_t(const tfbs_options_t& options,
                    const data_tfbs_t& data,
                    boost::optional<const alignment_set_t<>&> alignment_set =
                    boost::optional<const alignment_set_t<>&>());
         dpm_tfbs_t(const dpm_tfbs_t& dpm);
        ~dpm_tfbs_t();

        virtual dpm_tfbs_t* clone() const;

        friend void swap(dpm_tfbs_t& first, dpm_tfbs_t& second);

        // auxiliary types
        ////////////////////////////////////////////////////////////////////////
        typedef mixture_state_t::const_iterator cm_iterator;
        typedef cluster_t::const_iterator cl_iterator;
        typedef indexer_t::const_iterator da_iterator;

        // operators
        ////////////////////////////////////////////////////////////////////////
        friend std::ostream& operator<<(std::ostream& o, const dpm_tfbs_t& dpm);

        virtual dpm_tfbs_t& operator=(const mixture_model_t& mixture_model);

        // access methods
        ////////////////////////////////////////////////////////////////////////
        const dpm_tfbs_state_t& state() const;
              dpm_tfbs_state_t& state();

        const data_tfbs_t& data() const;

        const alignment_set_t<>& alignment_set() const;

        const dpm_tfbs_prior_t& process_prior() const;

        // methods
        ////////////////////////////////////////////////////////////////////////
        size_t mixture_components() const;
        size_t baseline_components() const;
        double background_mixture_weight(const range_t& range, cluster_t& cluster);
        double foreground_mixture_weight(const range_t& range, cluster_t& cluster);
        // override pure functions
        void   mixture_weights(const range_t& range, double log_weights[], cluster_tag_t tags[]) {
                // set default temperature
                mixture_weights(range, log_weights, tags, 1.0, true, -std::numeric_limits<double>::infinity());
        }
        void   mixture_weights(const range_t& range, double log_weights[], cluster_tag_t tags[],
                               const double temp,
                               const bool include_background = true,
                               const double baseline = -std::numeric_limits<double>::infinity());
        void   mixture_weights(const std::vector<range_t>& range_set, double log_weights[], cluster_tag_t cluster_tags[],
                               const double temp = 1.0,
                               const bool include_background = true);

        double likelihood() const;
        double posterior() const;

        // compute point estimates
        ////////////////////////////////////////////////////////////////////////
        dpm_partition_t map   (const sampling_history_t& histroy, bool optimize = false, bool verbose = false) const;
        dpm_partition_t mean  (const sampling_history_t& history, ssize_t take = -1, bool verbose = false) const;
        dpm_partition_t median(const sampling_history_t& history, ssize_t take = -1, bool verbose = false) const;

        // test methods
        ////////////////////////////////////////////////////////////////////////
        void test();
        void test_background();
        void test_moves();
        void test_posterior(cluster_t& cluster, const range_t& range);

protected:
        // baseline models
        std::vector<double> m_baseline_weights;
        std::vector<baseline_tag_t> m_baseline_tags;

        // data and clusters
        const data_tfbs_t* m_data;
        const alignment_set_t<>* m_alignment_set;
        dpm_tfbs_state_t m_state;

        // parameters
        double m_lambda;
        double m_lambda_log;
        double m_lambda_inv_log;

        // process priors
        dpm_tfbs_prior_t* m_process_prior;
};

#endif /* __TFBAYES_DPM_DPM_TFBS_HH__ */
