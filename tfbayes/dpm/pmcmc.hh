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

#ifndef PMCMC_HH
#define PMCMC_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <sampler.hh>

class population_mcmc_t : public sampler_t {
public:
        population_mcmc_t(size_t n);
        population_mcmc_t(const population_mcmc_t& pmcmc);
        virtual ~population_mcmc_t();

        population_mcmc_t* clone() const;

        // operators
        ////////////////////////////////////////////////////////////////////////
              sampler_t& operator[](size_t i)       { return *_population[i]; }
        const sampler_t& operator[](size_t i) const { return *_population[i]; }

        // methods
        ////////////////////////////////////////////////////////////////////////
        void sample(size_t n, size_t burnin);
        size_t size() const;

        const sampling_history_t& sampling_history() const;
        samples_t& samples();
        size_t sampling_steps() const;

protected:
        std::vector<sampler_t*> _population;
        const size_t _size;
        sampling_history_t* _sampling_history;
        samples_t* _samples;

        // private methods
        ////////////////////////////////////////////////////////////////////////
        void update_samples();
        void update_sampling_history();
};

#include <data-tfbs.hh>
#include <dpm-tfbs.hh>
#include <repl-server.hh>

class dpm_tfbs_sampler_t : public population_mcmc_t {
public:
        dpm_tfbs_sampler_t(const tfbs_options_t& options, const std::vector<std::string>& sequences, size_t n);
        virtual ~dpm_tfbs_sampler_t();

        data_tfbs_t _data;
        std::vector<dpm_tfbs_t*> _gdpm;

private:
        void _start_server();
        void _stop_server();

        std::string _socket_file;
        boost::asio::io_service _ios;
        server_t* _server;
        boost::thread* _bt;

        std::vector<save_queue_t<command_t*> > _command_queue;
        save_queue_t<std::string> _output_queue;
};

#endif /* PMCMC_HH */
