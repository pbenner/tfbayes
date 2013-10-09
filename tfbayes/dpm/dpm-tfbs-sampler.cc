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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/dpm/dpm-tfbs-sampler.hh>
#include <tfbayes/dpm/statistics.hh>

using namespace std;
using namespace boost::asio;
using namespace boost::asio::local;

// dpm_tfbs_sampler_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_sampler_t::dpm_tfbs_sampler_t(
        const tfbs_options_t& options,
        dpm_tfbs_t& dpm,
        dpm_tfbs_state_t& state,
        const indexer_t& indexer,
        const string name,
        save_queue_t<command_t*>& command_queue,
        save_queue_t<string>& output_queue,
        const sequence_data_t<data_tfbs_t::code_t>& phylogenetic_data)
        : gibbs_sampler_t(dpm, state, indexer, name),
          phylogenetic_data(&phylogenetic_data),
          _command_queue(&command_queue),
          _output_queue(&output_queue),
          _state(&state),
          _dpm(&dpm),
          _t0(options.initial_temperature)
{
        assert(options.initial_temperature >= 1.0);
}

dpm_tfbs_sampler_t::~dpm_tfbs_sampler_t()
{
        delete(_command_queue);
}

dpm_tfbs_sampler_t*
dpm_tfbs_sampler_t::clone() const {
        return new dpm_tfbs_sampler_t(*this);
}

// Gibbs samples
////////////////////////////////////////////////////////////////////////////////

bool
dpm_tfbs_sampler_t::_gibbs_sample(const index_i& index, const double temp, const bool optimize) {
        ////////////////////////////////////////////////////////////////////////
        // check if we can sample this element
        if (!_dpm->valid_for_sampling(index)) {
                return false;
        }
        ////////////////////////////////////////////////////////////////////////
        // release the element from its cluster
        cluster_tag_t old_cluster_tag = (*_state)[index];
        _state->remove(index, old_cluster_tag);
        size_t components = _dpm->mixture_components() + _dpm->baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        _dpm->mixture_weights(index, log_weights, cluster_tags, temp);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        cluster_tag_t new_cluster_tag;
        if (optimize) {
                new_cluster_tag = cluster_tags[select_max_component(components, log_weights)];
                if (new_cluster_tag != old_cluster_tag) {
                        flockfile(stdout);
                        cout << _name << ": "
                             << "moving " << (const seq_index_t&)index
                             << " from cluster " << old_cluster_tag
                             << " to cluster "   << new_cluster_tag
                             << endl;
                        fflush(stdout);
                        funlockfile(stdout);
                }
        }
        else {
                new_cluster_tag = cluster_tags[select_component(components, log_weights)];
        }

        ////////////////////////////////////////////////////////////////////////
        _state->add(index, new_cluster_tag);

        return old_cluster_tag != new_cluster_tag;
}

size_t
dpm_tfbs_sampler_t::_gibbs_sample(const double temp, const bool optimize) {
        size_t sum = 0;
        // the indexer needs to be constant since it is shared between
        // processes, so to shuffle the indices we first need to
        // obtain a copy
        vector<index_i*> indices(_indexer->sampling_begin(), _indexer->sampling_end());
        random_shuffle(indices.begin(), indices.end());
        // now sample
        for (vector<index_i*>::iterator it = indices.begin();
             it != indices.end(); it++) {
                if(_gibbs_sample(**it, temp, optimize)) sum+=1;
        }
        return sum;
}

// Gibbs block samples
////////////////////////////////////////////////////////////////////////////////

void
dpm_tfbs_sampler_t::_block_sample(cluster_t& cluster, const double temp, const bool optimize)
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

                _state->remove(range.index(), old_cluster_tag);
        }
        ////////////////////////////////////////////////////////////////////////
        // obtian the mixture probabilities
        size_t components = _dpm->mixture_components() + _dpm->baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        _dpm->mixture_weights(range_set, log_weights, cluster_tags, temp, false);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        cluster_tag_t new_cluster_tag;
        if (optimize) {
                new_cluster_tag = cluster_tags[select_max_component(components, log_weights)];
        }
        else {
                new_cluster_tag = cluster_tags[select_component(components, log_weights)];
        }

        ////////////////////////////////////////////////////////////////////////
        // print some information to stdout
        flockfile(stdout);
        if (_state[new_cluster_tag].size() != 0) {
                cout << _name << ": "
                     << "cluster " << old_cluster_tag
                     << " merged into cluster " << new_cluster_tag
                     << " (" << _state[new_cluster_tag].size()
                     << "+" << range_set.size() << ")." << endl;
        }
        fflush(stdout);
        funlockfile(stdout);

        ////////////////////////////////////////////////////////////////////////
        // move all elements to the new cluster
        for (vector<range_t>::iterator it = range_set.begin(); it != range_set.end(); it++)
        {
                const range_t& range = *it;

                _state->add(range.index(), new_cluster_tag);
        }
}

void
dpm_tfbs_sampler_t::_block_sample(const double temp, const bool optimize)
{
        // since clusters are modified it is not possible to simply
        // loop through the list of clusters, we need to be a bit more
        // careful here!
        vector<cluster_tag_t> used_clusters;
        for (cm_iterator it = _state->begin(); it != _state->end(); it++) {
                cluster_t& cluster = **it;
                if (cluster.cluster_tag() != _state->bg_cluster_tag) {
                        used_clusters.push_back(cluster.cluster_tag());
                }
        }

        // go through the list of used clusters and if they are still
        // used then generate a block sample
        for (vector<cluster_tag_t>::const_iterator it = used_clusters.begin(); it != used_clusters.end(); it++) {
                cluster_t& cluster = (*_state)[*it];
                if (cluster.size() != 0) {
                        _block_sample(cluster, temp, optimize);
                }
        }
}

// Metropolis-Hastings samples
////////////////////////////////////////////////////////////////////////////////

bool
dpm_tfbs_sampler_t::_metropolis_sample(cluster_t& cluster, const double temp) {
        double posterior_ref = _dpm->posterior();
        double posterior_tmp;
        stringstream ss;
        size_t size = cluster.size();

        if (_state->proposal(cluster, ss)) {
                posterior_tmp = _dpm->posterior();

                const double r = (double)rand()/RAND_MAX;
                /* posterior value is on log scale! */
                if (r <= min(exp((posterior_tmp - posterior_ref)/temp), 1.0)) {
                        goto accepted;
                }
                _state->restore();
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
dpm_tfbs_sampler_t::_metropolis_sample(const double temp) {
        for (cl_iterator it = _state->begin(); it != _state->end(); it++) {
                _metropolis_sample(**it, temp);
        }

        return true;
}

// Main
////////////////////////////////////////////////////////////////////////////////

size_t
dpm_tfbs_sampler_t::_sample(size_t i, size_t n, bool is_burnin) {
        // temperature for simulated annealing
        double temperature = 1.0;
        if (is_burnin) {
                // geometric decline of the temperature
                temperature = _t0*pow((1.0/_t0), (double)i/n);
        }
        flockfile(stdout);
        cout << _name << ": "
             << "temperature is " << temperature << endl;
        fflush(stdout);
        funlockfile(stdout);
        // save temperature
        _sampling_history.temperature[0].push_back(temperature);
        // call the standard hybrid sampler that first produces a
        // Gibbs sample and afterwards make a Metropolis-Hastings step
        size_t s = _gibbs_sample(temperature, false);
        _metropolis_sample(temperature);
        // do a Gibbs block sampling step, i.e. go through all
        // clusters and try to merge them
        _block_sample(temperature, false);
        // we are done with sampling here, now process commands
        flockfile(stdout);
        cout << _name << ": "
             << "Processing commands."
             << endl;
        fflush(stdout);
        funlockfile(stdout);
        while (!_command_queue->empty()) {
                stringstream ss;
                command_t* command = _command_queue->front();
                ss << command->operator()(*_state, *this);
                _command_queue->pop();
                delete(command);
                _output_queue->push(ss.str());
        }
        return s;
}

// Optimization
////////////////////////////////////////////////////////////////////////////////

bool
dpm_tfbs_sampler_t::optimize(cluster_t& cluster) {
        double posterior_ref = _dpm->posterior();
        double posterior_left;
        double posterior_right;
        stringstream ss;
        size_t size = cluster.size();

        _state->save(cluster.cluster_tag());
        _state->move_left(cluster);
        posterior_left = _dpm->posterior();
        _state->restore();

        _state->save(cluster.cluster_tag());
        _state->move_right(cluster);
        posterior_right = _dpm->posterior();
        _state->restore();

        if (posterior_left > posterior_ref && posterior_left > posterior_right) {
                _state->move_left(cluster);
                ss << "moved to the left";
        }
        else if (posterior_right > posterior_ref && posterior_right > posterior_left) {
                _state->move_right(cluster);
                ss << "moved to the right";
        }
        else {
                return false;
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

        return true;
}

void
dpm_tfbs_sampler_t::optimize() {
        // print old posterior value
        flockfile(stdout);
        cout << _name << ": "
             << "old posterior value: " << _dpm->posterior()
             << endl;
        fflush(stdout);
        funlockfile(stdout);
        // block optimization
        for (cl_iterator it = _state->begin(); it != _state->end(); it++) {
                optimize(**it);
        }
        // output resulting posterior value
        flockfile(stdout);
        cout << _name << ": "
             << "new posterior value: " << _dpm->posterior()
             << endl;
        fflush(stdout);
        funlockfile(stdout);
}

// dpm_tfbs_pmcmc_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_pmcmc_t::dpm_tfbs_pmcmc_t(
        const tfbs_options_t& options,
        const sampling_history_t& history)
        : population_mcmc_t(options.population_size),
          _options(options),
          _data(options.phylogenetic_file, options.tfbs_length),
          _alignment_set(options.alignment_file),
          _gdpm(options.population_size, NULL),
          _socket_file(options.socket_file),
          _server(NULL),
          _bt(NULL)
{
        const dpm_tfbs_t dpm_tfbs(options, _data, _alignment_set);

        for (size_t i = 0; i < _size; i++) {
                std::stringstream ss; ss << "Sampler " << i+1;
                save_queue_t<command_t*>* command_queue = new save_queue_t<command_t*>();
                _command_queue.push_back(command_queue);
                // clone dpm_tfbs instead of constructing a new one
                // for every sampler; this prevents that any
                // preprocessing is done several times
                _gdpm[i]       = dpm_tfbs.clone();
                _population[i] = new dpm_tfbs_sampler_t(options,
                                                        *_gdpm[i],
                                                        _gdpm[i]->state(),
                                                        _data,
                                                        ss.str(),
                                                        *command_queue,
                                                        _output_queue,
                                                        _data);
        }
        _start_server();
}

dpm_tfbs_pmcmc_t::~dpm_tfbs_pmcmc_t() {
        _stop_server();
}

const tfbs_options_t&
dpm_tfbs_pmcmc_t::options() const {
        return _options;
}

const data_tfbs_t&
dpm_tfbs_pmcmc_t::data() const {
        return _data;
}

const std::vector<dpm_tfbs_t*>&
dpm_tfbs_pmcmc_t::gdpm() const {
        return _gdpm;
}

dpm_partition_t
dpm_tfbs_pmcmc_t::map() const {
        return dpm_partition_t();
}

dpm_partition_t
dpm_tfbs_pmcmc_t::mean() const {
        return dpm_partition_t();
}

dpm_partition_t
dpm_tfbs_pmcmc_t::median() const {
        return dpm_partition_t();
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

void
dpm_tfbs_pmcmc_t::save(const string& filename) const
{
        if (filename == "") {
                cout << *this;
        }
        else {
                ofstream file;
                file.open(filename.c_str());
                file << *this;
                file.close();
        }
}

static
ostream& operator<< (ostream& o, const matrix<double>& m)
{
        for (size_t i = 0; i < m.size(); i++) {
                o << "\t";
                for (size_t j = 0; j < m[i].size(); j++) {
                        o << m[i][j] << " ";
                }
                o << endl;
        }
        return o;
}

ostream& operator<< (ostream& o, const dpm_tfbs_pmcmc_t& pmcmc)
{
        const sampling_history_t& history = pmcmc.sampling_history();

        o << "[Result]" << endl;

        o << "components =" << endl
          << std::noshowpoint
          << history.components;
        o << "switches =" << endl
          << std::noshowpoint
          << history.switches;
        o << "likelihood =" << endl
          << std::showpoint
          << history.likelihood;
        o << "posterior =" << endl
          << std::showpoint
          << history.posterior;
        o << "temperature =" << endl
          << std::showpoint
          << history.temperature;
        o << "partitions =" << endl
          << history.partitions;

        return o;
}
