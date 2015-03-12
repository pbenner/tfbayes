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

#include <iostream>
#include <sstream>
#include <cassert>

#include <boost/thread.hpp>

#include <tfbayes/dpm/pmcmc.hh>

using namespace std;

population_mcmc_t::population_mcmc_t(size_t n, const sampling_history_t& history)
        : m_population       (n, NULL)
        , m_size             (n)
        , m_sampling_history (history)
{
        assert(n >= 1);

        if (m_sampling_history.switches.size() == 0) {
                m_sampling_history.switches = matrix<double>(n, 0);
        }
        if (m_sampling_history.likelihood.size() == 0) {
                m_sampling_history.likelihood = matrix<double>(n, 0);
        }
        if (m_sampling_history.posterior.size() == 0) {
                m_sampling_history.posterior = matrix<double>(n, 0);
        }
        if (m_sampling_history.components.size() == 0) {
                m_sampling_history.components = matrix<double>(n, 0);
        }
        if (m_sampling_history.temperature.size() == 0) {
                m_sampling_history.temperature = matrix<double>(n, 0);
        }

        assert(m_sampling_history.switches   .size() == n);
        assert(m_sampling_history.likelihood .size() == n);
        assert(m_sampling_history.posterior  .size() == n);
        assert(m_sampling_history.components .size() == n);
        assert(m_sampling_history.temperature.size() == n);
}

population_mcmc_t::population_mcmc_t(const population_mcmc_t& pmcmc)
        : sampler_t    (pmcmc)
        , m_population (pmcmc.m_size, NULL)
        , m_size       (pmcmc.m_size)
{
        for (size_t i = 0; i < m_size; i++) {
                m_population[i] = (sampler_t*)pmcmc.m_population[i]->clone();
        }
}

population_mcmc_t::~population_mcmc_t()
{
        for (size_t i = 0; i < m_size; i++) {
                delete(m_population[i]);
        }
}

population_mcmc_t*
population_mcmc_t::clone() const {
        return new population_mcmc_t(*this);
}

population_mcmc_t&
population_mcmc_t::operator=(const sampler_t& sampler)
{
        population_mcmc_t tmp(static_cast<const population_mcmc_t&>(sampler));
        swap(*this, tmp);
        return *this;
}

void
swap(population_mcmc_t& first, population_mcmc_t& second) {
        swap(static_cast<sampler_t&>(first), static_cast<sampler_t&>(second));
        swap(first.m_population,       second.m_population);
        swap(first.m_size,             second.m_size);
        swap(first.m_sampling_history, second.m_sampling_history);
}

void
population_mcmc_t::update_sampling_history()
{
        // get sampling histories from population
        for (size_t i = 0; i < m_size; i++) {
                m_sampling_history.switches[i].insert(
                        m_sampling_history.switches[i].end(),
                        m_population[i]->sampling_history().switches[0].begin(),
                        m_population[i]->sampling_history().switches[0].end());
                m_sampling_history.likelihood[i].insert(
                        m_sampling_history.likelihood[i].end(),
                        m_population[i]->sampling_history().likelihood[0].begin(),
                        m_population[i]->sampling_history().likelihood[0].end());
                m_sampling_history.posterior[i].insert(
                        m_sampling_history.posterior[i].end(),
                        m_population[i]->sampling_history().posterior[0].begin(),
                        m_population[i]->sampling_history().posterior[0].end());
                m_sampling_history.components[i].insert(
                        m_sampling_history.components[i].end(),
                        m_population[i]->sampling_history().components[0].begin(),
                        m_population[i]->sampling_history().components[0].end());
                m_sampling_history.temperature[i].insert(
                        m_sampling_history.temperature[i].end(),
                        m_population[i]->sampling_history().temperature[0].begin(),
                        m_population[i]->sampling_history().temperature[0].end());
        }
        // copy cluster sizes
        for (size_t j = 0; j < m_population[0]->sampling_history().partitions.size(); j++) {
                for (size_t i = 0; i < m_size; i++) {
                        assert(j < m_population[i]->sampling_history().cluster_sizes.size());
                        m_sampling_history.cluster_sizes.push_back(
                                m_population[i]->sampling_history().cluster_sizes[j]);
                }
        }
        // copy partitions
        for (size_t j = 0; j < m_population[0]->sampling_history().partitions.size(); j++) {
                for (size_t i = 0; i < m_size; i++) {
                        assert(j < m_population[i]->sampling_history().partitions.size());
                        m_sampling_history.partitions.push_back(
                                m_population[i]->sampling_history().partitions[j]);
                }
        }
        // reset sampling history
        for (size_t i = 0; i < m_size; i++) {
                m_population[i]->sampling_history() = sampling_history_t();
        }
}

void
population_mcmc_t::operator()(size_t n, size_t burnin)
{
        std::vector<boost::thread*> threads(m_size);

        // sample
        for (size_t i = 0; i < m_size; i++) {
                threads[i] = new boost::thread(boost::ref(*m_population[i]), n, burnin);
        }
        // join threads
        for (size_t i = 0; i < m_size; i++) {
                threads[i]->join();
        }
        for (size_t i = 0; i < m_size; i++) {
                delete(threads[i]);
        }
        update_sampling_history();
}

const sampling_history_t&
population_mcmc_t::sampling_history() const {
        return m_sampling_history;
}

sampling_history_t&
population_mcmc_t::sampling_history() {
        return m_sampling_history;
}

size_t
population_mcmc_t::size() const {
        return m_size;
}
