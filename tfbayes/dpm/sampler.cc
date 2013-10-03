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

#include <sstream>

#include <gsl/gsl_randist.h>

#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/dpm/state.hh>
#include <tfbayes/dpm/statistics.hh>

using namespace std;

// Gibbs Sampler
////////////////////////////////////////////////////////////////////////////////

gibbs_sampler_t::gibbs_sampler_t(mixture_model_t& dpm,
                           gibbs_state_t& state,
                           const indexer_t& indexer,
                           const string name)
        : _dpm(&dpm),
          _name(name),
          _state(&state),
          _indexer(&indexer)
{
        // for sampling statistics
        _sampling_history.switches.   push_back(vector<double>());
        _sampling_history.likelihood. push_back(vector<double>());
        _sampling_history.posterior.  push_back(vector<double>());
        _sampling_history.components. push_back(vector<double>());
        _sampling_history.temperature.push_back(vector<double>());
}

gibbs_sampler_t::gibbs_sampler_t(const gibbs_sampler_t& sampler)
        : _dpm(sampler._dpm),
          _name(sampler._name),
          _state(sampler._state),
          _indexer(sampler._indexer),
          _sampling_history(*new sampling_history_t(sampler._sampling_history))
{
}

gibbs_sampler_t::~gibbs_sampler_t()
{
        delete(&_dpm);
        delete(&_sampling_history);
}

gibbs_sampler_t*
gibbs_sampler_t::clone() const {
        return new gibbs_sampler_t(*this);
}

bool
gibbs_sampler_t::_gibbs_sample(const index_i& index) {
        ////////////////////////////////////////////////////////////////////////
        // check if we can sample this element
        if (!_dpm.valid_for_sampling(index)) {
                return false;
        }
        ////////////////////////////////////////////////////////////////////////
        // release the element from its cluster
        cluster_tag_t old_cluster_tag = _state[index];
        _state.remove(index, old_cluster_tag);
        size_t components = _dpm.mixture_components() + _dpm.baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        _dpm.mixture_weights(index, log_weights, cluster_tags);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        cluster_tag_t new_cluster_tag = cluster_tags[select_component(components, log_weights)];

        ////////////////////////////////////////////////////////////////////////
        _state.add(index, new_cluster_tag);

        return old_cluster_tag != new_cluster_tag;
}

size_t
gibbs_sampler_t::_gibbs_sample() {
        size_t sum = 0;
        // the indexer needs to be constant since it is shared between
        // processes, so to shuffle the indices we first need to
        // obtain a copy
        vector<index_i*> indices(_indexer.sampling_begin(), _indexer.sampling_end());
        random_shuffle(indices.begin(), indices.end());
        // now sample
        for (vector<index_i*>::iterator it = indices.begin();
             it != indices.end(); it++) {
                if(_gibbs_sample(**it)) sum+=1;
        }
        return sum;
}

size_t
gibbs_sampler_t::_sample(size_t i, size_t n, bool is_burnin) {
        return _gibbs_sample();
}

const sampling_history_t&
gibbs_sampler_t::sampling_history() const {
        return _sampling_history;
}

void
gibbs_sampler_t::_update_sampling_history(size_t switches)
{
        _sampling_history.switches  [0].push_back(switches);
        _sampling_history.likelihood[0].push_back(_dpm.likelihood());
        _sampling_history.posterior [0].push_back(_dpm.posterior());
        _sampling_history.components[0].push_back(_dpm.mixture_components());
        _sampling_history.partitions   .push_back(_dpm.partition());
}

void
gibbs_sampler_t::sample(size_t n, size_t burnin) {
        // temperature for simulated annealing
        // burn in sampling
        for (size_t i = 0; i < burnin; i++) {
                flockfile(stdout);
                cout << _name << ": "
                     << "Burn in... [" << i+1 << "]"
                     << "[ Cluster: " << _state << "]"
                     << endl;
                fflush(stdout);
                funlockfile(stdout);
                _update_sampling_history(_sample(i, burnin, true));
        }
        // sample `n' times
        for (size_t i = 0; i < n; i++) {
                // loop through all elements
                flockfile(stdout);
                cout << _name << ": "
                     << "Sampling... [" << i+1 << "]"
                     << "[ Cluster: " << _state << "]"
                     << endl;
                fflush(stdout);
                funlockfile(stdout);
                _update_sampling_history(_sample(i, n, false));
        }
}

string
gibbs_sampler_t::name() const {
        return _name;
}
