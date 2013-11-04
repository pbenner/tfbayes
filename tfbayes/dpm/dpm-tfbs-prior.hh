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

#ifndef __TFBAYES_DPM_DPM_TFBS_PRIOR_HH__
#define __TFBAYES_DPM_DPM_TFBS_PRIOR_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/dpm/cluster.hh>
#include <tfbayes/dpm/dpm-tfbs-state.hh>
#include <tfbayes/utility/clonable.hh>

class dpm_tfbs_prior_t : public clonable {
public:
        virtual ~dpm_tfbs_prior_t() {}

        dpm_tfbs_prior_t* clone() const = 0;

        virtual double log_predictive(const cluster_t& cluster, const dpm_tfbs_state_t& state) const = 0;
        virtual double joint(const dpm_tfbs_state_t& state) const = 0;
};

class pitman_yor_prior : public dpm_tfbs_prior_t {
public:
        pitman_yor_prior(double alpha, double discount,
                         cluster_tag_t bg_cluster_tag);

        pitman_yor_prior* clone() const;

        double log_predictive(const cluster_t& cluster, const dpm_tfbs_state_t& state) const;
        double joint(const dpm_tfbs_state_t& state) const;

protected:
        typedef mixture_state_t::const_iterator cl_iterator;

        const double alpha;
        const double discount;
        const cluster_tag_t bg_cluster_tag;
};

class uniform_prior : public dpm_tfbs_prior_t {
public:
        uniform_prior(double alpha);

        uniform_prior* clone() const;

        double log_predictive(const cluster_t& cluster, const dpm_tfbs_state_t& state) const;
        double joint(const dpm_tfbs_state_t& state) const;

protected:
        const double alpha;
};

class poppe_prior : public dpm_tfbs_prior_t {
public:
        poppe_prior();

        poppe_prior* clone() const;

        double log_predictive(const cluster_t& cluster, const dpm_tfbs_state_t& state) const;
        double joint(const dpm_tfbs_state_t& state) const;
};

#endif /* __TFBAYES_DPM_DPM_TFBS_PRIOR_HH__ */
