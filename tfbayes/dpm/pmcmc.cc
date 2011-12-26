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

#include <iostream>
#include <sstream>

#include <assert.h>
#include <pthread.h>
#include <unistd.h>

#include <pmcmc.hh>

using namespace std;

PopulationMCMC::PopulationMCMC(Sampler& sampler, size_t n)
        : _population(n, NULL), _size(n),
          _sampling_history(NULL), _samples(NULL)
{
        assert(n >= 1);

        _population[0] = &sampler;
        _population[0]->set_name("Sampler 1");
        for (size_t i = 1; i < _size; i++) {
                stringstream ss;
                ss << "Sampler " << i+1;
                _population[i] = (Sampler*)sampler.clone();
                _population[i]->set_name(ss.str());
        }
}

PopulationMCMC::PopulationMCMC(const PopulationMCMC& pmcmc)
        : _population(pmcmc._size, NULL), _size(pmcmc._size),
          _sampling_history(NULL), _samples(NULL)
{
        for (size_t i = 0; i < _size; i++) {
                _population[i] = (Sampler*)pmcmc._population[i]->clone();
        }
}

PopulationMCMC::~PopulationMCMC()
{
        for (size_t i = 0; i < _size; i++) {
                delete(_population[i]);
        }
        if (_sampling_history != NULL) {
                delete(_sampling_history);
        }
        if (_samples != NULL) {
                delete(_samples);
        }
}

PopulationMCMC*
PopulationMCMC::clone() const {
        return new PopulationMCMC(*this);
}

void
PopulationMCMC::update_sampling_history()
{
        if (_sampling_history != NULL) {
                delete(_sampling_history);
        }
        _sampling_history = new sampling_history_t();
        for (size_t i = 0; i < _size; i++) {
                _sampling_history->switches.push_back(_population[i]->sampling_history().switches[0]);
                _sampling_history->likelihood.push_back(_population[i]->sampling_history().likelihood[0]);
                _sampling_history->components.push_back(_population[i]->sampling_history().components[0]);
        }
}

void
PopulationMCMC::update_samples()
{
        if (_samples != NULL) {
                delete(_samples);
        }
        // allocate memory
        _samples = new samples_t();
        for (size_t i = 0; i < _population[0]->samples().probabilities.size(); i++) {
                _samples->probabilities.push_back(vector<double>(_population[0]->samples().probabilities[i].size(), 0));
        }
        // loop through samples
        for (size_t i = 0; i < _population[0]->samples().probabilities.size(); i++) {
                for (size_t j = 0; j < _population[0]->samples().probabilities[i].size(); j++) {
                        // average over population
                        double sum = 0;
                        for (size_t k = 0; k < _size; k++) {
                                sum += _population[k]->samples().probabilities[i][j];
                        }
                        // save average
                        _samples->probabilities[i][j] = sum/(double)_size;
                }
        }
        // graph updates
        for (size_t i = 0; i < _size; i++) {
                _samples->graph.insert(_population[i]->samples().graph);
                _population[i]->samples().graph.cleanup();
        }
}

typedef struct {
        Sampler* sampler;
        size_t n;
        size_t burnin;
} pthread_data_t;

static
void * sample_thread(void* _data)
{
        pthread_data_t* data  = (pthread_data_t*)_data;
        Sampler* sampler      = data->sampler;
        const size_t n        = data->n;
        const size_t burnin   = data->burnin;

        sampler->sample(n, burnin);

        return NULL;
}

void
PopulationMCMC::sample(size_t n, size_t burnin)
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
        update_samples();
}

const sampling_history_t&
PopulationMCMC::sampling_history() const {
        return *_sampling_history;
}

samples_t&
PopulationMCMC::samples() {
        return *_samples;
}

size_t
PopulationMCMC::size() const {
        return _size;
}

size_t
PopulationMCMC::sampling_steps() const {
        return _population[0]->sampling_steps();
}

void
PopulationMCMC::set_name(const string name)
{
}
