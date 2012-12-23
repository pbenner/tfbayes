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

#ifndef DPM_TFBS_COMMAND_HH
#define DPM_TFBS_COMMAND_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <string>

class dpm_tfbs_sampler_t;
class dpm_tfbs_state_t;

class command_t {
public:
        virtual ~command_t() {};

        virtual std::string operator()(const dpm_tfbs_state_t& state, dpm_tfbs_sampler_t& sampler) const = 0;
};

class print_cluster_counts_t : public command_t {
public:
        print_cluster_counts_t(ssize_t cluster_tag);

        std::string operator()(const dpm_tfbs_state_t& state, dpm_tfbs_sampler_t& sampler) const;

private:
        ssize_t _cluster_tag;
};

class print_cluster_elements_t : public command_t {
public:
        print_cluster_elements_t(ssize_t cluster_tag);

        std::string operator()(const dpm_tfbs_state_t& state, dpm_tfbs_sampler_t& sampler) const;

private:
        ssize_t _cluster_tag;
};

class print_likelihood_t : public command_t {
public:
        std::string operator()(const dpm_tfbs_state_t& state, dpm_tfbs_sampler_t& sampler) const;

};

class print_posterior_t : public command_t {
public:
        std::string operator()(const dpm_tfbs_state_t& state, dpm_tfbs_sampler_t& sampler) const;

};

#endif /* DPM_TFBS_COMMAND_HH */
