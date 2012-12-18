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

#include <sstream>

#include <gsl/gsl_randist.h>

#include <sampler.hh>
#include <state.hh>
#include <statistics.hh>

using namespace std;

// Gibbs Sampler
////////////////////////////////////////////////////////////////////////////////

gibbs_sampler_t::gibbs_sampler_t(mixture_model_t& dpm,
                           gibbs_state_t& state,
                           const indexer_t& indexer,
                           const string name)
        : _dpm(dpm),
          _name(name),
          _state(state),
          _indexer(indexer),
          _sampling_steps(0),
          _sampling_history(*new sampling_history_t())
{
        // for sampling statistics
        _sampling_history.switches.push_back(vector<double>());
        _sampling_history.likelihood.push_back(vector<double>());
        _sampling_history.posterior.push_back(vector<double>());
        _sampling_history.components.push_back(vector<size_t>());
        _sampling_history.switches[0].push_back(0);
        _sampling_history.likelihood[0].push_back(_dpm.likelihood());
        _sampling_history.posterior[0].push_back(_dpm.posterior());
        _sampling_history.components[0].push_back(1);
}

gibbs_sampler_t::gibbs_sampler_t(const gibbs_sampler_t& sampler)
        : _dpm(sampler._dpm),
          _name(sampler._name),
          _state(sampler._state),
          _indexer(sampler._indexer),
          _sampling_steps(0),
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

size_t
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
        for (indexer_t::sampling_iterator it = _indexer.sampling_begin();
             it != _indexer.sampling_end(); it++) {
                if(_gibbs_sample(**it)) sum+=1;
        }
        return sum;
}

bool
gibbs_sampler_t::_sample() {
        return _gibbs_sample();
}

const sampling_history_t&
gibbs_sampler_t::sampling_history() const {
        return _sampling_history;
}

samples_t&
gibbs_sampler_t::samples() {
        return _dpm.samples();
}

size_t
gibbs_sampler_t::sampling_steps() const {
        return _sampling_steps;
}

void
gibbs_sampler_t::sample(size_t n, size_t burnin) {
        double sum;
        double switches;
        // burn in sampling
        for (size_t i = 0; i < burnin; i++) {
                flockfile(stdout);
                cout << _name << ": "
                     << "Burn in... [" << i+1 << "]"
                     << "[ Cluster: " << _state << "]"
                     << endl;
                fflush(stdout);
                funlockfile(stdout);
                sum      = _sample();
                switches = _indexer.elements() > 0 ? sum/(double)_indexer.elements() : 0;
                _sampling_history.likelihood[0].push_back(_dpm.likelihood());
                _sampling_history.posterior[0].push_back(_dpm.posterior());
                _sampling_history.components[0].push_back(_dpm.mixture_components());
                _sampling_history.switches[0].push_back(switches);
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
                sum      = _sample();
                switches = _indexer.elements() > 0 ? sum/(double)_indexer.elements() : 0;
                _sampling_history.likelihood[0].push_back(_dpm.likelihood());
                _sampling_history.posterior[0].push_back(_dpm.posterior());
                _sampling_history.components[0].push_back(_dpm.mixture_components());
                _sampling_history.switches[0].push_back(switches);
                _dpm.update_samples(_sampling_steps);
                _sampling_steps++;
        }
}

string
gibbs_sampler_t::name() const {
        return _name;
}

// Hybrid Sampler (combines a Gibbs sampler with a Metropolis-Hastings sampler)
////////////////////////////////////////////////////////////////////////////////

hybrid_sampler_t::hybrid_sampler_t(
        mixture_model_t& dpm,
        hybrid_state_t& state,
        const indexer_t& indexer,
        const string name,
        bool optimize)
        : gibbs_sampler_t(dpm, state, indexer, name),
          _state(state), _optimize(optimize)
{
}

hybrid_sampler_t*
hybrid_sampler_t::clone() const {
        return new hybrid_sampler_t(*this);
}

bool
hybrid_sampler_t::_sample() {
        size_t s = _gibbs_sample();
                   _metropolis_sample();
        return s;
}

bool
hybrid_sampler_t::_metropolis_sample(cluster_t& cluster) {
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
hybrid_sampler_t::_metropolis_sample() {
        for (cl_iterator it = _state.begin(); it != _state.end(); it++) {
                _metropolis_sample(**it);
        }

        return true;
}
