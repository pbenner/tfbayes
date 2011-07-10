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

#include <datatypes.hh>
#include <dpm.hh>

class GibbsSampler {
public:
         GibbsSampler(DPM& dpm, const Data& data);
        ~GibbsSampler();

        const sampling_history_t& sampling_history() const;

        void sample(size_t n, size_t burnin);
private:
        // private methods
        bool _sample(const element_t& element);

        // the mixture model
        DPM& _dpm;
        const Data& _data;

        // gibbs sampler history
        size_t _sampling_steps;
        sampling_history_t& _sampling_history;
};

#endif /* SAMPLER_HH */
