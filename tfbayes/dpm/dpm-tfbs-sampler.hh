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

#ifndef DPM_TFBS_SAMPLER_HH
#define DPM_TFBS_SAMPLER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/dpm/dpm-tfbs-state.hh>
#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/dpm/utility.hh>

class command_t;

class dpm_tfbs_sampler_t : public gibbs_sampler_t {
public:
        dpm_tfbs_sampler_t(
                dpm_tfbs_t& dpm,
                dpm_tfbs_state_t& state,
                const indexer_t& indexer,
                const std::string name,
                bool optimize,
                save_queue_t<command_t*>& command_queue,
                save_queue_t<std::string>& output_queue,
                const sequence_data_t<data_tfbs_t::code_t>& sequences);

        // after sampling this function locally optimizes the current state
        void optimize(cluster_t& cluster);
        void optimize();

        dpm_tfbs_sampler_t* clone() const;

        // auxiliary types
        ////////////////////////////////////////////////////////////////////////
        typedef mixture_state_t::const_iterator cm_iterator;
        typedef mixture_state_t::iterator cl_iterator;

        const sequence_data_t<data_tfbs_t::code_t> sequences;

protected:
        bool _sample();
        void _block_sample();
        void _block_sample(cluster_t& cluster);
        bool _metropolis_sample();
        bool _metropolis_sample(cluster_t& cluster);

        save_queue_t<command_t*>& _command_queue;
        save_queue_t<std::string>& _output_queue;
        dpm_tfbs_state_t& _state;
        dpm_tfbs_t& _dpm;
        const bool _optimize;
};

#include <pmcmc.hh>
#include <dpm-tfbs.hh>
#include <dpm-tfbs-command.hh>
#include <dpm-tfbs-repl.hh>

class dpm_tfbs_pmcmc_t : public population_mcmc_t {
public:
        dpm_tfbs_pmcmc_t(
                const tfbs_options_t& options,
                const sequence_data_t<data_tfbs_t::code_t>& sequences,
                size_t n);
        virtual ~dpm_tfbs_pmcmc_t();

        data_tfbs_t _data;
        std::vector<dpm_tfbs_t*> _gdpm;

        void optimize();

private:
        void _start_server();
        void _stop_server();

        std::string _socket_file;
        boost::asio::io_service _ios;
        server_t* _server;
        boost::thread* _bt;

        std::vector<save_queue_t<command_t*>* > _command_queue;
        save_queue_t<std::string> _output_queue;
        const sequence_data_t<data_tfbs_t::code_t> _sequences;
};

#endif /* DPM_TFBS_SAMPLER_HH */
