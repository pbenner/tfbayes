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

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <tfbayes/dpm/dpm-tfbs-sampler.hh>
#include <tfbayes/utility/statistics.hh>

using namespace std;
using namespace boost::asio;
using namespace boost::asio::local;

// dpm_tfbs_sampler_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_sampler_t::dpm_tfbs_sampler_t(
        const tfbs_options_t& options,
        const dpm_tfbs_t& dpm_tfbs,
        const data_tfbs_t& data,
        save_queue_t<string>& output_queue)
        :
        gibbs_sampler_t  (dpm_tfbs, data),
        phylogenetic_data(&data),
        _output_queue    (&output_queue),
        _t0              (options.initial_temperature)
{
        assert(options.initial_temperature >= 1.0);
}

dpm_tfbs_sampler_t::dpm_tfbs_sampler_t(const dpm_tfbs_sampler_t& sampler)
        : gibbs_sampler_t  (sampler),
          phylogenetic_data(sampler.phylogenetic_data),
          _output_queue    (sampler._output_queue),
          _t0              (sampler._t0)
{ }

dpm_tfbs_sampler_t::~dpm_tfbs_sampler_t()
{
}

void
swap(dpm_tfbs_sampler_t& first, dpm_tfbs_sampler_t& second) {
        swap(static_cast<gibbs_sampler_t&>(first), static_cast<gibbs_sampler_t&>(second));
        swap(first.phylogenetic_data, second.phylogenetic_data);
        swap(first._command_queue,    second._command_queue);
        swap(first._output_queue,     second._output_queue);
        swap(first._t0,               second._t0);
}

dpm_tfbs_sampler_t*
dpm_tfbs_sampler_t::clone() const {
        return new dpm_tfbs_sampler_t(*this);
}

dpm_tfbs_sampler_t&
dpm_tfbs_sampler_t::operator=(const sampler_t& sampler)
{
        dpm_tfbs_sampler_t tmp(static_cast<const dpm_tfbs_sampler_t&>(sampler));
        swap(*this, tmp);
        return *this;
}

save_queue_t<command_t*>&
dpm_tfbs_sampler_t::command_queue()
{
        return _command_queue;
}

const save_queue_t<command_t*>&
dpm_tfbs_sampler_t::command_queue() const
{
        return _command_queue;
}

dpm_tfbs_t&
dpm_tfbs_sampler_t::dpm()
{
        return *static_cast<dpm_tfbs_t*>(_dpm);
}

const dpm_tfbs_t&
dpm_tfbs_sampler_t::dpm() const
{
        return *static_cast<const dpm_tfbs_t*>(_dpm);
}

// Gibbs samples
////////////////////////////////////////////////////////////////////////////////

bool
dpm_tfbs_sampler_t::_gibbs_sample(const index_i& index, const double temp) {
        ////////////////////////////////////////////////////////////////////////
        // check if we can sample this element
        if (!dpm().state().valid_tfbs_position(index)) {
                return false;
        }
        ////////////////////////////////////////////////////////////////////////
        // first release the element from its cluster
        cluster_tag_t old_cluster_tag = dpm().state()[index];
        state().remove(index, old_cluster_tag);
        ////////////////////////////////////////////////////////////////////////
        size_t components = dpm().mixture_components() + dpm().baseline_components();
        double log_weights1[components];
        double log_weights2[components];
        cluster_tag_t cluster_tags1[components];
        cluster_tag_t cluster_tags2[components];
        range_t range1(index, dpm().state().tfbs_length, false);
        range_t range2(index, dpm().state().tfbs_length, true );
        ////////////////////////////////////////////////////////////////////////
        // compute weights
        dpm().mixture_weights(range1, log_weights1, cluster_tags1, temp);
        dpm().mixture_weights(range2, log_weights2, cluster_tags2, temp, log_weights1[components-1]);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        std::pair<size_t, size_t> result = select_component2(components, log_weights1, log_weights2, gen());
        ////////////////////////////////////////////////////////////////////////
        if (result.first == 1) {
                dpm().state().add(range1, cluster_tags1[result.second]);
                return old_cluster_tag != cluster_tags1[result.second];
        }
        else {
                dpm().state().add(range2, cluster_tags2[result.second]);
                return old_cluster_tag != cluster_tags2[result.second];
        }
}

size_t
dpm_tfbs_sampler_t::_gibbs_sample(const double temp) {
        size_t sum = 0;
        // the indexer needs to be constant since it is shared between
        // processes, so to shuffle the indices we first need to
        // obtain a copy
        vector<index_i*> indices(_indexer->sampling_begin(), _indexer->sampling_end());
        random_shuffle(indices.begin(), indices.end());
        // now sample
        for (vector<index_i*>::iterator it = indices.begin();
             it != indices.end(); it++) {
                if(_gibbs_sample(**it, temp)) sum+=1;
        }
        return sum;
}

// Gibbs block samples
////////////////////////////////////////////////////////////////////////////////

void
dpm_tfbs_sampler_t::_block_sample(cluster_t& cluster, const double temp)
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

                dpm().state().remove(range, old_cluster_tag);
        }
        ////////////////////////////////////////////////////////////////////////
        // obtian the mixture probabilities
        size_t components = dpm().mixture_components() + dpm().baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        dpm().mixture_weights(range_set, log_weights, cluster_tags, temp, false);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        cluster_tag_t new_cluster_tag = cluster_tags[select_component(components, log_weights, gen())];

        ////////////////////////////////////////////////////////////////////////
        // print some information to stdout
        flockfile(stdout);
        if (state()[new_cluster_tag].size() != 0) {
                cout << boost::format("%s: cluster %d merged with cluster %d (%d + %d)")
                        % _name % old_cluster_tag % new_cluster_tag
                        % state()[new_cluster_tag].size() % range_set.size()
                     << endl;
        }
        fflush(stdout);
        funlockfile(stdout);

        ////////////////////////////////////////////////////////////////////////
        // move all elements to the new cluster
        for (size_t i = 0; i < range_set.size(); i++) {
                dpm().state().add(range_set[i], new_cluster_tag);
        }
}

void
dpm_tfbs_sampler_t::_block_sample(const double temp)
{
        // since clusters are modified it is not possible to simply
        // loop through the list of clusters, we need to be a bit more
        // careful here!
        vector<cluster_tag_t> used_clusters;
        for (cm_iterator it = dpm().state().begin(); it != dpm().state().end(); it++) {
                cluster_t& cluster = **it;
                if (cluster.cluster_tag() != dpm().state().bg_cluster_tag) {
                        used_clusters.push_back(cluster.cluster_tag());
                }
        }

        // go through the list of used clusters and if they are still
        // used then generate a block sample
        for (vector<cluster_tag_t>::const_iterator it = used_clusters.begin(); it != used_clusters.end(); it++) {
                cluster_t& cluster = state()[*it];
                if (cluster.size() != 0) {
                        _block_sample(cluster, temp);
                }
        }
}

// Metropolis-Hastings samples
////////////////////////////////////////////////////////////////////////////////

bool
dpm_tfbs_sampler_t::_metropolis_sample(cluster_t& cluster, const double temp) {
        double posterior_ref = dpm().posterior();
        double posterior_tmp;
        stringstream ss;
        size_t size = cluster.size();

        /* allocate a uniform distribution on the unit inverval */
        boost::random::uniform_01<> dist;

        if (dpm().state().proposal(cluster, ss, gen())) {
                posterior_tmp = dpm().posterior();

                /* posterior value is on log scale! */
                if (dist(gen()) <= min(exp((posterior_tmp - posterior_ref)/temp), 1.0)) {
                        goto accepted;
                }
                dpm().state().restore();
        }
        return false;

accepted:
        flockfile(stdout);
        cout << boost::format("%s: cluster %d: %s accepted (%d -> %d)")
                % _name % cluster.cluster_tag() % ss.str() % size % cluster.size()
             << endl;
        fflush(stdout);
        funlockfile(stdout);

        return true;
}

bool
dpm_tfbs_sampler_t::_metropolis_sample(const double temp) {
        for (cl_iterator it = dpm().state().begin(); it != dpm().state().end(); it++) {
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
        size_t s = _gibbs_sample(temperature);
        _metropolis_sample(temperature);
        // do a Gibbs block sampling step, i.e. go through all
        // clusters and try to merge them
        _block_sample(temperature);
        // we are done with sampling here, now process commands
        flockfile(stdout);
        cout << boost::format("%s: Processing commands.") % _name
             << endl;
        fflush(stdout);
        funlockfile(stdout);
        while (!_command_queue.empty()) {
                stringstream ss;
                command_t* command = _command_queue.front();
                ss << command->operator()(dpm().state(), *this);
                _command_queue.pop();
                delete(command);
                _output_queue->push(ss.str());
        }
        return s;
}

// dpm_tfbs_pmcmc_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_pmcmc_t::dpm_tfbs_pmcmc_t(
        const tfbs_options_t& options,
        const sampling_history_t& history)
        : population_mcmc_t(options.population_size, history),
          _options(options),
          _data(options.phylogenetic_file, options.tfbs_length),
          _alignment_set(options.alignment_file),
          _socket_file(options.socket_file),
          _server(NULL),
          _bt(NULL)
{
        for (size_t i = 0; i < _size; i++) {
                // initialize dpm_tfbs
                dpm_tfbs_t dpm_tfbs(options, _data, _alignment_set);
                // check if we have to resume from an old state
                if (history.partitions.size() != 0) {
                        if (history.partitions.size() % _size != 0) {
                                cerr << "Cannot resume from previous sampling run: Number of"
                                     << " partitions does not match."
                                     << endl;
                                exit(EXIT_FAILURE);
                        }
                        // resume partition
                        const dpm_partition_t& partition = history.partitions
                                [history.partitions.size()-_size+i];
                        dpm_tfbs.state().set_partition(partition);
                }
                // initialize sampler
                _population[i] = new dpm_tfbs_sampler_t(options, dpm_tfbs, _data, _output_queue);
                // make some noise
                std::stringstream ss;
                ss << "Sampler " << i+1;
                operator[](i).name() = ss.str();
        }
        _start_server();
}

dpm_tfbs_pmcmc_t::dpm_tfbs_pmcmc_t(const dpm_tfbs_pmcmc_t& sampler)
        : population_mcmc_t(sampler),
          _data(sampler._data)
{
        cerr << "Cannot copy dpm_tfbs_pmcmc_t!"
             << endl;
        exit(EXIT_FAILURE);
}

dpm_tfbs_pmcmc_t::~dpm_tfbs_pmcmc_t() {
        _stop_server();
}

dpm_tfbs_pmcmc_t*
dpm_tfbs_pmcmc_t::clone() const {
        return new dpm_tfbs_pmcmc_t(*this);
}

void
swap(dpm_tfbs_pmcmc_t& first, dpm_tfbs_pmcmc_t& second) {
        swap(static_cast<population_mcmc_t&>(first),
             static_cast<population_mcmc_t&>(second));
        swap(first._options,       second._options);
        swap(first._data,          second._data);
        swap(first._alignment_set, second._alignment_set);
        swap(first._socket_file,   second._socket_file);
        swap(first._server,        second._server);
        swap(first._bt,            second._bt);
        swap(first._output_queue,  second._output_queue);
}

dpm_tfbs_pmcmc_t&
dpm_tfbs_pmcmc_t::operator=(const sampler_t& sampler)
{
        dpm_tfbs_pmcmc_t tmp(static_cast<const dpm_tfbs_pmcmc_t&>(sampler));
        swap(*this, tmp);
        return *this;
}

const tfbs_options_t&
dpm_tfbs_pmcmc_t::options() const {
        return _options;
}

const data_tfbs_t&
dpm_tfbs_pmcmc_t::data() const {
        return _data;
}

void
dpm_tfbs_pmcmc_t::_start_server() {
        vector<save_queue_t<command_t*>* > command_queue;
        for (size_t i = 0; i < _size; i++) {
                command_queue.push_back(&operator[](i).command_queue());
        }
        if (_socket_file != "" && _server == NULL) {
                remove(_socket_file.c_str());
                _server = new server_t(_ios, _socket_file, command_queue, _output_queue);
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
