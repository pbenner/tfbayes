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

#ifndef SAMPLER_HH
#define SAMPLER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <clonable.hh>
#include <datatypes.hh>
#include <indexer.hh>
#include <mixture-model.hh>
#include <state.hh>

class Sampler : public clonable {
public:
        virtual void sample(size_t n, size_t burnin) = 0;
        virtual const sampling_history_t& sampling_history() const = 0;
        virtual samples_t& samples() = 0;
        virtual size_t sampling_steps() const = 0;
};

class GibbsSampler : public Sampler {
public:
         GibbsSampler(mixture_model_t& dpm,
                      gibbs_state_t& state,
                      const Indexer& indexer);
         GibbsSampler(const GibbsSampler& sampler);
        ~GibbsSampler();

        GibbsSampler* clone() const;

        void sample(size_t n, size_t burnin);
        const sampling_history_t& sampling_history() const;
        samples_t& samples();
        size_t sampling_steps() const;

protected:
        // private methods
        virtual bool _sample();
        size_t _gibbs_sample(const index_i& index);
        size_t _gibbs_sample();
        // the mixture model
        mixture_model_t& _dpm;

private:
        gibbs_state_t& _state;
        const Indexer& _indexer;

        // gibbs sampler history
        size_t _sampling_steps;
        sampling_history_t& _sampling_history;
};

class HybridSampler : public GibbsSampler {
public:
        HybridSampler(mixture_model_t& dpm,
                      hybrid_state_t& state,
                      const Indexer& indexer);

protected:
        virtual bool _sample();
                bool _metropolis_sample();

private:
        typedef mixture_state_t::iterator cl_iterator;

        hybrid_state_t& _state;
};

#endif /* SAMPLER_HH */
