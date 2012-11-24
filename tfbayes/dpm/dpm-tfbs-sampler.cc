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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <dpm-tfbs-sampler.hh>

using namespace std;
using namespace boost::asio;
using namespace boost::asio::local;

// dpm_tfbs_sampler_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_sampler_t::dpm_tfbs_sampler_t(
        mixture_model_t& dpm,
        dpm_tfbs_state_t& state,
        const indexer_t& indexer,
        const string name,
        bool optimize,
        save_queue_t<command_t*>& command_queue,
        save_queue_t<string>& output_queue,
        const sequence_data_t<data_tfbs_t::code_t>& sequences)
        : hybrid_sampler_t(dpm, state, indexer, name, optimize),
          sequences(sequences),
          _command_queue(command_queue),
          _output_queue(output_queue),
          _state(state) {
}

bool
dpm_tfbs_sampler_t::_sample() {
        size_t s = hybrid_sampler_t::_sample();
        flockfile(stdout);
        cout << _name << ": "
             << "Processing commands."
             << endl;
        fflush(stdout);
        funlockfile(stdout);
        while (!_command_queue.empty()) {
                stringstream ss;
                command_t* command = _command_queue.front();
                ss << command->operator()(_state, *this);
                _command_queue.pop();
                delete(command);
                _output_queue.push(ss.str());
        }
        return s;
}

// dpm_tfbs_pmcmc_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_pmcmc_t::dpm_tfbs_pmcmc_t(
        const tfbs_options_t& options,
        const sequence_data_t<data_tfbs_t::code_t>& sequences,
        size_t n)
        : population_mcmc_t(n),
          _data(sequences, options.tfbs_length),
          _gdpm(n, NULL),
          _socket_file(options.socket_file),
          _server(NULL),
          _bt(NULL),
          _sequences(sequences) {
        for (size_t i = 0; i < _size; i++) {
                std::stringstream ss; ss << "Sampler " << i+1;
                save_queue_t<command_t*>* command_queue = new save_queue_t<command_t*>();
                _command_queue.push_back(command_queue);
                _gdpm[i]       = new dpm_tfbs_t(options, _data);
                _population[i] = new dpm_tfbs_sampler_t(*_gdpm[i], _gdpm[i]->state(),
                                                        _data,
                                                        ss.str(),
                                                        options.metropolis_optimize,
                                                        *command_queue,
                                                        _output_queue,
                                                        _sequences);
        }
        update_samples();
        _start_server();
}

dpm_tfbs_pmcmc_t::~dpm_tfbs_pmcmc_t() {
        _stop_server();
}

void
dpm_tfbs_pmcmc_t::_start_server() {
        if (_socket_file != "" && _server == NULL) {
                remove(_socket_file.c_str());
                _server = new server_t(_ios, _socket_file, _command_queue, _output_queue);
                _bt     = new boost::thread(boost::bind(&io_service::run, &_ios));
        }
}

void
dpm_tfbs_pmcmc_t::_stop_server() {
        if (_socket_file != "" && _server != NULL) {
                _ios.stop();
                _bt->join();
                delete(_bt);
                delete(_server);
                remove(_socket_file.c_str());
        }
}
