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

#ifndef __TFBAYES_DPM_DPM_TFBS_SAMPLER_HH__
#define __TFBAYES_DPM_DPM_TFBS_SAMPLER_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <sstream>

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/dpm/dpm-tfbs-state.hh>
#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/dpm/save-queue.hh>

class command_t;

class dpm_tfbs_sampler_t : public gibbs_sampler_t {
public:
        dpm_tfbs_sampler_t(
                const tfbs_options_t& options,
                const dpm_tfbs_t& dpm,
                const data_tfbs_t& data,
                save_queue_t<std::string>& output_queue
                );
         dpm_tfbs_sampler_t(const dpm_tfbs_sampler_t& sampler);

        ~dpm_tfbs_sampler_t();

        dpm_tfbs_sampler_t* clone() const;

        friend void swap(dpm_tfbs_sampler_t& first, dpm_tfbs_sampler_t& second);

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual dpm_tfbs_sampler_t& operator=(const sampler_t& sampler);

        // after sampling this function locally optimizes the current state
        ////////////////////////////////////////////////////////////////////////
        bool optimize(cluster_t& cluster);
        void optimize();

        // access methods
        ////////////////////////////////////////////////////////////////////////
        const save_queue_t<command_t*>& command_queue() const;
              save_queue_t<command_t*>& command_queue();
        const dpm_tfbs_t& dpm() const;
              dpm_tfbs_t& dpm();

        // auxiliary types
        ////////////////////////////////////////////////////////////////////////
        typedef mixture_state_t::const_iterator cm_iterator;
        typedef mixture_state_t::iterator cl_iterator;

        // phylogenetic and alignment data for diagnostics
        ////////////////////////////////////////////////////////////////////////
        const sequence_data_t<data_tfbs_t::code_t>* phylogenetic_data;

protected:
        size_t m_sample(size_t i, size_t n, bool is_burnin);
        size_t m_sample(size_t i, size_t n, double temp, bool optimize);
        using gibbs_sampler_t::m_gibbs_sample;
        bool   m_gibbs_sample(const index_t& index, double temp, bool optimize);
        size_t m_gibbs_sample(double temp = 1.0, bool optimize = false);
        void m_block_sample(double temp, bool optimize);
        void m_block_sample(cluster_tag_t cluster_tag, double temp, bool optimize);
        bool m_metropolis_proposal_size(cluster_t& cluster, std::stringstream& ss);
        bool m_metropolis_proposal_move(cluster_t& cluster, std::stringstream& ss);
        bool m_metropolis_sample(double temp, bool optimize);
        bool m_metropolis_sample(cluster_tag_t cluster_tag, double temp, bool optimize,
                                 boost::function<bool (cluster_t& cluster, std::stringstream& ss)> f);
        void m_update_sampling_history(size_t switches);
        save_queue_t<command_t*> m_command_queue;
        save_queue_t<std::string>* m_output_queue;

        // initial temperature for simulated annealing
        double m_t0;
        bool   m_block_samples;
        size_t m_block_samples_period;
        size_t m_metropolis_proposals;
        bool   m_optimize;
        size_t m_optimize_period;
        size_t m_verbose;
};

#include <pmcmc.hh>
#include <dpm-tfbs.hh>
#include <dpm-tfbs-command.hh>
#include <dpm-tfbs-repl.hh>

class dpm_tfbs_pmcmc_t : public population_mcmc_t {
public:
        dpm_tfbs_pmcmc_t(
                const tfbs_options_t& options,
                const sampling_history_t& history = sampling_history_t());
        dpm_tfbs_pmcmc_t(const dpm_tfbs_pmcmc_t& sampler);

        virtual ~dpm_tfbs_pmcmc_t();

        dpm_tfbs_pmcmc_t* clone() const;

        friend void swap(dpm_tfbs_pmcmc_t& first, dpm_tfbs_pmcmc_t& second);

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual       dpm_tfbs_sampler_t& operator[](size_t i)
                { return *static_cast<      dpm_tfbs_sampler_t*>(m_population[i]); }
        virtual const dpm_tfbs_sampler_t& operator[](size_t i) const
                { return *static_cast<const dpm_tfbs_sampler_t*>(m_population[i]); }

        dpm_tfbs_pmcmc_t& operator=(const sampler_t& sampler);

        // access methods
        ////////////////////////////////////////////////////////////////////////
        const tfbs_options_t& options() const;
        const data_tfbs_t& data() const;

        // save sampling results
        ////////////////////////////////////////////////////////////////////////
        void save(const std::string& filename) const;

protected:
        void m_start_server();
        void m_stop_server();

        tfbs_options_t m_options;

        data_tfbs_t m_data;
        alignment_set_t<> m_alignment_set;

        std::string m_socket_file;
        boost::asio::io_service m_ios;
        server_t* m_server;
        boost::thread* m_bt;

        save_queue_t<std::string> m_output_queue;
};

std::ostream& operator<< (std::ostream& o, const dpm_tfbs_pmcmc_t& pmcmc);

#endif /* __TFBAYES_DPM_DPM_TFBS_SAMPLER_HH__ */
