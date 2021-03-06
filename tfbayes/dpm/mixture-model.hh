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

#ifndef __TFBAYES_DPM_MIXTURE_MODEL_HH__
#define __TFBAYES_DPM_MIXTURE_MODEL_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <mixture-state.hh>
#include <dpm-partition.hh>

#include <tfbayes/utility/clonable.hh>

class mixture_model_t : public virtual clonable {
public:
        virtual ~mixture_model_t() {}

        mixture_model_t* clone() const = 0;

        // methods
        ////////////////////////////////////////////////////////////////////////
        virtual size_t mixture_components() const = 0;
        virtual size_t baseline_components() const = 0;
        // compute cumulative mixture weights on log scale (not normalized!)
        virtual void   mixture_weights(const range_t& range, double log_weights[], cluster_tag_t tags[]) = 0;
        virtual double likelihood() const = 0;
        // compute the posterior value of the current cluster
        // assignment (usually not normalized)
        virtual double posterior() const = 0;
        virtual mixture_state_t& state() = 0;
        virtual const mixture_state_t& state() const = 0;
};

#endif /* __TFBAYES_DPM_MIXTURE_MODEL_HH__ */
