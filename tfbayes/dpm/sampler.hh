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

class sampler_t : public clonable {
public:
        virtual void sample(size_t n, size_t burnin) = 0;
        virtual const sampling_history_t& sampling_history() const = 0;
        virtual samples_t& samples() = 0;
        virtual size_t sampling_steps() const = 0;
        virtual std::string name() const { return std::string(); }
protected:
};

class gibbs_sampler_t : public sampler_t {
public:
         gibbs_sampler_t(mixture_model_t& dpm,
                         gibbs_state_t& state,
                         const indexer_t& indexer,
                         const std::string name = "");
         gibbs_sampler_t(const gibbs_sampler_t& sampler);
        ~gibbs_sampler_t();

        gibbs_sampler_t* clone() const;

        void sample(size_t n, size_t burnin);
        const sampling_history_t& sampling_history() const;
        samples_t& samples();
        size_t sampling_steps() const;
        std::string name() const;

protected:
        // private methods
        virtual bool _sample();
        size_t _gibbs_sample(const index_i& index);
        size_t _gibbs_sample();
        // the mixture model
        mixture_model_t& _dpm;
        std::string _name;

private:
        gibbs_state_t& _state;
        const indexer_t& _indexer;

        // gibbs sampler history
        size_t _sampling_steps;
        sampling_history_t& _sampling_history;
};

class hybrid_sampler_t : public gibbs_sampler_t {
public:
        hybrid_sampler_t(mixture_model_t& dpm,
                         hybrid_state_t& state,
                         const indexer_t& indexer,
                         const std::string name = "",
                         bool optimize = false);

        hybrid_sampler_t* clone() const;

protected:
        virtual bool _sample();
                bool _metropolis_sample();
                bool _metropolis_sample(cluster_t& cluster);
private:
        typedef mixture_state_t::iterator cl_iterator;

        hybrid_state_t& _state;
        const bool _optimize;
};

#endif /* SAMPLER_HH */
