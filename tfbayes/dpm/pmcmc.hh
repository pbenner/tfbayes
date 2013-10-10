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

#ifndef PMCMC_HH
#define PMCMC_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <tfbayes/dpm/sampling-history.hh>
#include <tfbayes/dpm/sampler.hh>

class population_mcmc_t : public sampler_t {
public:
        population_mcmc_t(size_t n, const sampling_history_t& history);
        population_mcmc_t(const population_mcmc_t& pmcmc);
        virtual ~population_mcmc_t();

        population_mcmc_t* clone() const;

        friend void swap(population_mcmc_t& first, population_mcmc_t& second);

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual       sampler_t& operator[](size_t i)       { return *_population[i]; }
        virtual const sampler_t& operator[](size_t i) const { return *_population[i]; }

        virtual population_mcmc_t& operator=(const sampler_t& sampler);

        // methods
        ////////////////////////////////////////////////////////////////////////
        void sample(size_t n, size_t burnin);
        size_t size() const;

        const sampling_history_t& sampling_history() const;

protected:
        std::vector<sampler_t*> _population;
        size_t _size;
        sampling_history_t _sampling_history;

        // protected methods
        ////////////////////////////////////////////////////////////////////////
        virtual void update_sampling_history();
};

#endif /* PMCMC_HH */
