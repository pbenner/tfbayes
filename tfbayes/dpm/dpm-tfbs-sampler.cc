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
#include <statistics.hh>

using namespace std;
using namespace boost::asio;
using namespace boost::asio::local;

// dpm_tfbs_sampler_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_sampler_t::dpm_tfbs_sampler_t(
        dpm_tfbs_t& dpm,
        dpm_tfbs_state_t& state,
        const indexer_t& indexer,
        const string name,
        bool optimize,
        save_queue_t<command_t*>& command_queue,
        save_queue_t<string>& output_queue,
        const sequence_data_t<data_tfbs_t::code_t>& sequences)
        : gibbs_sampler_t(dpm, state, indexer, name),
          sequences(sequences),
          _command_queue(command_queue),
          _output_queue(output_queue),
          _state(state),
          _dpm(dpm),
          _optimize(optimize)
{ }

dpm_tfbs_sampler_t*
dpm_tfbs_sampler_t::clone() const {
        return new dpm_tfbs_sampler_t(*this);
}

void
dpm_tfbs_sampler_t::_block_sample(cluster_t& cluster)
{
        vector<range_t> range_set;
        cluster_tag_t old_cluster_tag = cluster.cluster_tag();

        ////////////////////////////////////////////////////////////////////////
        // fill range_set
        for (cluster_t::iterator it = cluster.begin(); it != cluster.end(); it++)
        {
                const range_t& range = *it;

                range_set.push_back(range);
        }

        ////////////////////////////////////////////////////////////////////////
        // release all elemente from the cluster
        for (vector<range_t>::iterator it = range_set.begin(); it != range_set.end(); it++)
        {
                const range_t& range = *it;

                _state.remove(range.index(), old_cluster_tag);
        }
        ////////////////////////////////////////////////////////////////////////
        // obtian the mixture probabilities
        size_t components = _dpm.mixture_components() + _dpm.baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        _dpm.mixture_weights(range_set, log_weights, cluster_tags, false);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        cluster_tag_t new_cluster_tag = cluster_tags[select_component(components, log_weights)];

        ////////////////////////////////////////////////////////////////////////
        // print some information to stdout
        if (_state[new_cluster_tag].size() != 0) {
                cout << "Cluster " << old_cluster_tag
                     << " is merged into cluster " << new_cluster_tag
                     << " (" << _state[new_cluster_tag].size()
                     << "+" << range_set.size() << ")." << endl;
        }

        ////////////////////////////////////////////////////////////////////////
        // move all elements to the new cluster
        for (vector<range_t>::iterator it = range_set.begin(); it != range_set.end(); it++)
        {
                const range_t& range = *it;

                _state.add(range.index(), new_cluster_tag);
        }

}

void
dpm_tfbs_sampler_t::_block_sample()
{
        // since clusters are modified it is not possible to simply
        // loop through the list of clusters, we need to be a bit more
        // careful here!
        vector<cluster_tag_t> used_clusters;
        for (cm_iterator it = _state.begin(); it != _state.end(); it++) {
                cluster_t& cluster = **it;
                if (cluster.cluster_tag() != _state.bg_cluster_tag) {
                        used_clusters.push_back(cluster.cluster_tag());
                }
        }

        // go through the list of used clusters and if they are still
        // used then generate a block sample
        for (vector<cluster_tag_t>::const_iterator it = used_clusters.begin(); it != used_clusters.end(); it++) {
                cluster_t& cluster = _state[*it];
                if (cluster.size() != 0) {
                        _block_sample(cluster);
                }
        }
}

bool
dpm_tfbs_sampler_t::_metropolis_sample(cluster_t& cluster) {
        double posterior_ref = _dpm.posterior();
        double posterior_tmp;
        stringstream ss;
        size_t size = cluster.size();

        if (_state.proposal(cluster, ss)) {
                posterior_tmp = _dpm.posterior();

                if (_optimize && posterior_tmp > posterior_ref) {
                        goto accepted;
                }
                if (!_optimize) {
                        const double r = (double)rand()/RAND_MAX;
                        /* posterior value is on log scale! */
                        if (r <= min(exp(posterior_tmp - posterior_ref), 1.0)) {
                                goto accepted;
                        }
                }
                _state.restore();
        }
        return false;

accepted:
        flockfile(stdout);
        cout << _name << ": "
             << "cluster "
             << cluster.cluster_tag()
             << ": "
             << ss.str()
             << " accepted "
             << "("  << size
             << "->" << cluster.size()
             << ")"
             << endl;
        fflush(stdout);
        funlockfile(stdout);

        return true;
}

bool
dpm_tfbs_sampler_t::_metropolis_sample() {
        for (cl_iterator it = _state.begin(); it != _state.end(); it++) {
                _metropolis_sample(**it);
        }

        return true;
}

bool
dpm_tfbs_sampler_t::_sample() {
        // call the standard hybrid sampler that first produces a
        // Gibbs sample and afterwards make a Metropolis-Hastings step
        size_t s = _gibbs_sample();
                   _metropolis_sample();
        // do a Gibbs block sampling step, i.e. go through all
        // clusters and try to merge them
        _block_sample();
        // we are done with sampling here, now process commands
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

void
dpm_tfbs_sampler_t::optimize(cluster_t& cluster) {
        double posterior_ref = _dpm.posterior();
        double posterior_left;
        double posterior_right;
        stringstream ss;
        size_t size = cluster.size();

        _state.save(cluster.cluster_tag());
        _state.move_left(cluster);
        posterior_left = _dpm.posterior();
        _state.restore();

        _state.save(cluster.cluster_tag());
        _state.move_right(cluster);
        posterior_right = _dpm.posterior();
        _state.restore();

        if (posterior_left > posterior_ref && posterior_left > posterior_right) {
                _state.move_left(cluster);
                ss << "moved to the left";
        }
        else if (posterior_right > posterior_ref && posterior_right > posterior_left) {
                _state.move_right(cluster);
                ss << "moved to the right";
        }
        else {
                return;
        }

        flockfile(stdout);
        cout << _name << ": "
             << "cluster "
             << cluster.cluster_tag()
             << " "
             << ss.str()
             << " "
             << "("  << size
             << "->" << cluster.size()
             << ")"
             << endl;
        fflush(stdout);
        funlockfile(stdout);
}

void
dpm_tfbs_sampler_t::optimize() {
        for (cl_iterator it = _state.begin(); it != _state.end(); it++) {
                optimize(**it);
        }
        _dpm.update_samples(_sampling_steps);
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
dpm_tfbs_pmcmc_t::optimize() {
        for (size_t i = 0; i < _size; i++) {
                static_cast<dpm_tfbs_sampler_t *>(_population[i])->optimize();
        }
        update_samples();
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
