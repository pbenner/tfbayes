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

#include <pmcmc.hh>

using namespace std;

PopulationMCMC::PopulationMCMC(const Sampler& sampler, size_t n)
        : _population(n, NULL), _size(n),
          _sampling_history(*new sampling_history_t())
{
        for (size_t i = 0; i < _size; i++) {
                _population[i] = (Sampler*)sampler.clone();
        }
}

PopulationMCMC::PopulationMCMC(const PopulationMCMC& pmcmc)
        : _population(pmcmc._size, NULL), _size(pmcmc._size),
          _sampling_history(*new sampling_history_t())
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
        for (size_t i = 0; i < _size; i++) {
                _population[i]->sample(n, burnin);
        }
}

const sampling_history_t&
PopulationMCMC::sampling_history() const {
        return _sampling_history;
}

size_t
PopulationMCMC::size() const {
        return _size;
}
