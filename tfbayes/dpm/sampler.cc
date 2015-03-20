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
#include <sys/time.h>

#include <boost/thread/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <tfbayes/dpm/sampler.hh>
#include <tfbayes/dpm/state.hh>
#include <tfbayes/utility/statistics.hh>

using namespace std;

// Sampler
////////////////////////////////////////////////////////////////////////////////

sampler_t::sampler_t(const string& name)
        : m_name (name)
        , m_gen  () {
        /* sleep for a millisecond to make sure that we get
         * a unique seed */
        boost::this_thread::sleep(boost::posix_time::milliseconds(1));
        /* seed generator */
        struct timeval tv;
        gettimeofday(&tv, NULL);
        m_gen.seed(tv.tv_sec*tv.tv_usec);
}

// Gibbs Sampler
////////////////////////////////////////////////////////////////////////////////

gibbs_sampler_t::gibbs_sampler_t(const mixture_model_t& dpm,
                                 const indexer_t& indexer,
                                 const string name,
                                 bool verbose)
        : sampler_t  (name)
        , m_dpm     (dpm.clone())
        , m_indexer (&indexer)
        , m_verbose (verbose)
{
        // for sampling statistics
        m_sampling_history.switches.   push_back(vector<double>());
        m_sampling_history.likelihood. push_back(vector<double>());
        m_sampling_history.posterior.  push_back(vector<double>());
        m_sampling_history.components. push_back(vector<double>());
        m_sampling_history.temperature.push_back(vector<double>());
}

gibbs_sampler_t::gibbs_sampler_t(const gibbs_sampler_t& sampler)
        : sampler_t          (sampler)
        , m_dpm              (sampler.m_dpm->clone())
        , m_indexer          (sampler.m_indexer)
        , m_sampling_history (sampler.m_sampling_history)
        , m_verbose          (sampler.m_verbose)
{
}

gibbs_sampler_t::~gibbs_sampler_t()
{
        delete(m_dpm);
}

void
swap(gibbs_sampler_t& first, gibbs_sampler_t& second)
{
        swap(static_cast<sampler_t&>(first), static_cast<sampler_t&>(second));
        swap(first.m_dpm,              second.m_dpm);
        swap(first.m_name,             second.m_name);
        swap(first.m_indexer,          second.m_indexer);
        swap(first.m_sampling_history, second.m_sampling_history);
        swap(first.m_verbose,          second.m_verbose);
}

gibbs_sampler_t*
gibbs_sampler_t::clone() const {
        return new gibbs_sampler_t(*this);
}

gibbs_sampler_t&
gibbs_sampler_t::operator=(const sampler_t& sampler)
{
        gibbs_sampler_t tmp(static_cast<const gibbs_sampler_t&>(sampler));
        swap(*this, tmp);
        return *this;
}

bool
gibbs_sampler_t::m_gibbs_sample(const index_t& index)
{
        range_t range(index,1);
        ////////////////////////////////////////////////////////////////////////
        // release the element from its cluster
        cluster_tag_t old_cluster_tag = state()[index];
        state().remove(range);
        size_t components = m_dpm->mixture_components() + m_dpm->baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        m_dpm->mixture_weights(range, log_weights, cluster_tags);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        cluster_tag_t new_cluster_tag = cluster_tags[select_component(components, log_weights, gen())];

        ////////////////////////////////////////////////////////////////////////
        state().add(range, new_cluster_tag);

        return old_cluster_tag != new_cluster_tag;
}

size_t
gibbs_sampler_t::m_gibbs_sample() {
        size_t sum = 0;
        // the indexer needs to be constant since it is shared between
        // processes, so to shuffle the indices we first need to
        // obtain a copy
        vector<index_t> indices(m_indexer->sampling_begin(), m_indexer->sampling_end());
        random_shuffle(indices.begin(), indices.end());
        // now sample
        for (vector<index_t>::const_iterator it = indices.begin();
             it != indices.end(); it++) {
                if(m_gibbs_sample(*it)) sum+=1;
        }
        return sum;
}

size_t
gibbs_sampler_t::m_sample(size_t i, size_t n, bool is_burnin) {
        return m_gibbs_sample();
}

const sampling_history_t&
gibbs_sampler_t::sampling_history() const {
        return m_sampling_history;
}

sampling_history_t&
gibbs_sampler_t::sampling_history() {
        return m_sampling_history;
}

const mixture_model_t&
gibbs_sampler_t::dpm() const
{
        return *m_dpm;
}

mixture_model_t&
gibbs_sampler_t::dpm()
{
        return *m_dpm;
}

const gibbs_state_t&
gibbs_sampler_t::state() const {
        return static_cast<const gibbs_state_t&>(m_dpm->state());
}

gibbs_state_t&
gibbs_sampler_t::state() {
        return static_cast<gibbs_state_t&>(m_dpm->state());
}

void
gibbs_sampler_t::m_update_sampling_history(size_t switches)
{
        vector<double> tmp;

        BOOST_FOREACH(const cluster_t* cluster, m_dpm->state()) {
                tmp.push_back(cluster->size());
        }
        m_sampling_history.switches  [0].push_back(switches);
        m_sampling_history.likelihood[0].push_back(m_dpm->likelihood());
        m_sampling_history.posterior [0].push_back(m_dpm->posterior());
        m_sampling_history.components[0].push_back(m_dpm->mixture_components());
        m_sampling_history.partitions   .push_back(m_dpm->state().partition());
        m_sampling_history.cluster_sizes.push_back(tmp);
}

void
gibbs_sampler_t::operator()(size_t n, size_t burnin) {
        // temperature for simulated annealing
        // burn in sampling
        for (size_t i = 0; i < burnin; i++) {
                if (m_verbose) {
                        flockfile(stderr);
                        cerr << m_name << ": "
                             << "Burnin step " << i+1 << ":" << endl
                             << state() << endl;
                        fflush(stderr);
                        funlockfile(stderr);
                }
                m_update_sampling_history(m_sample(i, burnin, true));
        }
        // sample `n' times
        for (size_t i = 0; i < n; i++) {
                // loop through all elements
                if (m_verbose) {
                        flockfile(stderr);
                        cerr << m_name << ": "
                             << "Sampling step " << i+1 << ":" << endl
                             << state() << endl;
                        fflush(stderr);
                        funlockfile(stderr);
                }
                m_update_sampling_history(m_sample(i, n, false));
        }
}
