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

#include <assert.h>
#include <pthread.h>
#include <unistd.h>

#include <tfbayes/dpm/pmcmc.hh>

using namespace std;

population_mcmc_t::population_mcmc_t(size_t n)
        : _population(n, NULL),
          _size(n)
{
        assert(n >= 1);
}

population_mcmc_t::population_mcmc_t(const population_mcmc_t& pmcmc)
        : _population(pmcmc._size, NULL),
          _size(pmcmc._size)
{
        for (size_t i = 0; i < _size; i++) {
                _population[i] = (sampler_t*)pmcmc._population[i]->clone();
        }
}

population_mcmc_t::~population_mcmc_t()
{
        for (size_t i = 0; i < _size; i++) {
                delete(_population[i]);
        }
}

population_mcmc_t*
population_mcmc_t::clone() const {
        return new population_mcmc_t(*this);
}

void
population_mcmc_t::update_sampling_history()
{
        // reset sampling history
        _sampling_history = sampling_history_t();
        // get sampling histories from population
        for (size_t i = 0; i < _size; i++) {
                _sampling_history.switches   .push_back(_population[i]->sampling_history().switches   [0]);
                _sampling_history.likelihood .push_back(_population[i]->sampling_history().likelihood [0]);
                _sampling_history.posterior  .push_back(_population[i]->sampling_history().posterior  [0]);
                _sampling_history.components .push_back(_population[i]->sampling_history().components [0]);
                _sampling_history.temperature.push_back(_population[i]->sampling_history().temperature[0]);
        }
        for (size_t j = 0; j < _population[0]->sampling_history().partitions.size(); j++) {
                for (size_t i = 0; i < _size; i++) {
                        assert(j < _population[i]->sampling_history().partitions.size());
                        _sampling_history.partitions.push_back(
                                _population[i]->sampling_history().partitions[j]);
                }
        }
}

typedef struct {
        sampler_t* sampler;
        size_t n;
        size_t burnin;
} pthread_data_t;

static
void * sample_thread(void* _data)
{
        pthread_data_t* data  = (pthread_data_t*)_data;
        sampler_t* sampler    = data->sampler;
        const size_t n        = data->n;
        const size_t burnin   = data->burnin;

        sampler->sample(n, burnin);

        return NULL;
}

void
population_mcmc_t::sample(size_t n, size_t burnin)
{
        pthread_data_t data[_size];
        pthread_t threads[_size];
        int rc;

        for (size_t i = 0; i < _size; i++) {
                data[i].sampler = _population[i];
                data[i].n       = n;
                data[i].burnin  = burnin;
        }

        // sample
        for (size_t i = 0; i < _size; i++) {
                rc = pthread_create(&threads[i], NULL, sample_thread, (void *)&data[i]);
                if (rc) {
                        cerr << "Couldn't create thread." << endl;
                        exit(EXIT_FAILURE);
                }
        }
        // join threads
        for (size_t i = 0; i < _size; i++) {
                rc = pthread_join(threads[i], NULL);
                if (rc) {
                        cerr << "Couldn't join thread." << endl;
                        exit(EXIT_FAILURE);
                }
        }
        update_sampling_history();
}

const sampling_history_t&
population_mcmc_t::sampling_history() const {
        return _sampling_history;
}

size_t
population_mcmc_t::size() const {
        return _size;
}
