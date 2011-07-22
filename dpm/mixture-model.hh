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

#ifndef MIXTURE_MODEL_HH
#define MIXTURE_MODEL_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <clustermanager.hh>

class DPM : public clonable {
public:
        virtual ~DPM() {}

        virtual DPM* clone() const = 0;

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual       Cluster& operator[](cluster_tag_t c)       = 0;
        virtual const Cluster& operator[](cluster_tag_t c) const = 0;

        // methods
        ////////////////////////////////////////////////////////////////////////
        virtual size_t mixture_components() const = 0;
        virtual void   mixture_weights(const index_t& index, double weights[], cluster_tag_t tags[]) = 0;
        virtual void   add(const index_t& index, cluster_tag_t tag) = 0;
        virtual void   remove(const index_t& index, cluster_tag_t tag) = 0;
        virtual void   update_posterior(size_t sampling_steps) = 0;
        virtual double likelihood() const = 0;
        virtual bool   valid_for_sampling(const index_t& index) const = 0;
        virtual const posterior_t& posterior() const = 0;
        virtual const ClusterManager& cluster_manager() const = 0;
};

#endif /* MIXTURE_MODEL_HH */
