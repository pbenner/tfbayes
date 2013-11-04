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

#ifndef __TFBAYES_DPM_SAMPLER_HH__
#define __TFBAYES_DPM_SAMPLER_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/random/mersenne_twister.hpp>

#include <tfbayes/dpm/datatypes.hh>
#include <tfbayes/dpm/dpm-sampling-history.hh>
#include <tfbayes/dpm/indexer.hh>
#include <tfbayes/dpm/mixture-model.hh>
#include <tfbayes/dpm/state.hh>
#include <tfbayes/utility/clonable.hh>

class sampler_t : public clonable {
public:
        sampler_t(const std::string& name = "");
        sampler_t(const sampler_t& sampler)
                : _name(sampler._name),
                  _gen (sampler._gen) { }

        friend void swap(sampler_t& first, sampler_t& second) {
                std::swap(first._name, second._name);
                std::swap(first._gen,  second._gen);
        }

        virtual sampler_t& operator=(const sampler_t& sampler) = 0;

        virtual void operator()(size_t n, size_t burnin) = 0;

        virtual const sampling_history_t& sampling_history() const = 0;
        virtual       sampling_history_t& sampling_history() = 0;
        virtual const std::string& name() const { return _name; }
        virtual       std::string& name()       { return _name; }
        virtual boost::random::mt19937& gen()   { return _gen;  }
protected:
        std::string _name;
        boost::random::mt19937 _gen;
};

class gibbs_sampler_t : public sampler_t {
public:
         gibbs_sampler_t(const mixture_model_t& dpm,
                         const indexer_t& indexer,
                         const std::string name = "");
         gibbs_sampler_t(const gibbs_sampler_t& sampler);
        ~gibbs_sampler_t();

        friend void swap(gibbs_sampler_t& first, gibbs_sampler_t& second);

        gibbs_sampler_t* clone() const;

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual gibbs_sampler_t& operator=(const sampler_t& sampler);

        // methods
        ////////////////////////////////////////////////////////////////////////
        void operator()(size_t n, size_t burnin);

        // access methods
        ////////////////////////////////////////////////////////////////////////
        virtual const sampling_history_t& sampling_history() const;
        virtual       sampling_history_t& sampling_history();
        virtual const mixture_model_t& dpm() const;
        virtual       mixture_model_t& dpm();
        virtual const gibbs_state_t& state() const;
        virtual       gibbs_state_t& state();

protected:
        // private methods
        virtual size_t _sample(size_t i, size_t n, bool is_burnin);
        virtual bool   _gibbs_sample(const index_i& index);
        virtual size_t _gibbs_sample();
        virtual void   _update_sampling_history(size_t switches);
        // the mixture model
        mixture_model_t* _dpm;

        const indexer_t* _indexer;

        // gibbs sampler history
        sampling_history_t _sampling_history;
};

#endif /* __TFBAYES_DPM_SAMPLER_HH__ */
