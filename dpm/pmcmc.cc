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

#include <assert.h>

#include <pmcmc.hh>

using namespace std;

PopulationMCMC::PopulationMCMC(Sampler* sampler, size_t n)
        : _population(n, NULL), _size(n),
          _sampling_history(NULL)
{
        assert(n >= 1);

        _population[0] = sampler;
        for (size_t i = 1; i < _size; i++) {
                _population[i] = (Sampler*)sampler->clone();
        }
}

PopulationMCMC::PopulationMCMC(const PopulationMCMC& pmcmc)
        : _population(pmcmc._size, NULL), _size(pmcmc._size),
          _sampling_history(NULL)
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
}

PopulationMCMC*
PopulationMCMC::clone() const {
        return new PopulationMCMC(*this);
}

void
PopulationMCMC::sample(size_t n, size_t burnin)
{
        // sample
        for (size_t i = 0; i < _size; i++) {
                _population[i]->sample(n, burnin);
        }

        // update sampling history
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

const sampling_history_t&
PopulationMCMC::sampling_history() const {
        return *_sampling_history;
}

const Model&
PopulationMCMC::model() const {
        return _model;
}

const posterior_t&
PopulationMCMC::posterior() const {
        return _population[0]->posterior();
}

size_t
PopulationMCMC::size() const {
        return _size;
}
