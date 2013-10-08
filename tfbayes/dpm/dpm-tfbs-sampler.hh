/* Copyright (C) 2011, 2012 Philipp Benner
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

#ifndef DPM_TFBS_SAMPLER_HH
#define DPM_TFBS_SAMPLER_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/dpm/dpm-tfbs-state.hh>
#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/dpm/utility.hh>

class command_t;

class dpm_tfbs_sampler_t : public gibbs_sampler_t {
public:
        dpm_tfbs_sampler_t(
                const tfbs_options_t& options,
                dpm_tfbs_t& dpm,
                dpm_tfbs_state_t& state,
                const indexer_t& indexer,
                const std::string name,
                save_queue_t<command_t*>& command_queue,
                save_queue_t<std::string>& output_queue,
                const sequence_data_t<data_tfbs_t::code_t>& phylogenetic_data);

        // after sampling this function locally optimizes the current state
        bool optimize(cluster_t& cluster);
        void optimize();

        dpm_tfbs_sampler_t* clone() const;

        // phylogenetic and alignment data for diagnostics
        ////////////////////////////////////////////////////////////////////////
        const sequence_data_t<data_tfbs_t::code_t> phylogenetic_data;

        // auxiliary types
        ////////////////////////////////////////////////////////////////////////
        typedef mixture_state_t::const_iterator cm_iterator;
        typedef mixture_state_t::iterator cl_iterator;

protected:
        size_t _sample(size_t i, size_t n, bool is_burnin);
        using gibbs_sampler_t::_gibbs_sample;
        bool   _gibbs_sample(const index_i& index, const double temp, const bool optimize);
        size_t _gibbs_sample(const double temp = 1.0, const bool optimize = false);
        void _block_sample(const double temp, const bool optimize = false);
        void _block_sample(cluster_t& cluster, const double temp, const bool optimize);
        bool _metropolis_sample(const double temp);
        bool _metropolis_sample(cluster_t& cluster, const double temp);

        save_queue_t<command_t*>& _command_queue;
        save_queue_t<std::string>& _output_queue;
        dpm_tfbs_state_t& _state;
        dpm_tfbs_t& _dpm;

        // initial temperature for simulated annealing
        const double _t0;
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
        virtual ~dpm_tfbs_pmcmc_t();

        const tfbs_options_t& options() const;
        const data_tfbs_t& data() const;
        const std::vector<dpm_tfbs_t*>& gdpm() const;

        // compute point estimates
        dpm_partition_t map() const;
        dpm_partition_t mean() const;
        dpm_partition_t median() const;
        void save(const std::string& filename) const;

protected:
        void _start_server();
        void _stop_server();

        const tfbs_options_t* _options;

        data_tfbs_t _data;
        alignment_set_t<short> _alignment_set;
        std::vector<dpm_tfbs_t*> _gdpm;

        std::string _socket_file;
        boost::asio::io_service _ios;
        server_t* _server;
        boost::thread* _bt;

        std::vector<save_queue_t<command_t*>* > _command_queue;
        save_queue_t<std::string> _output_queue;
};

std::ostream& operator<< (std::ostream& o, const dpm_tfbs_pmcmc_t& pmcmc);

#endif /* DPM_TFBS_SAMPLER_HH */
