/* Copyright (C) 2011-2015 Philipp Benner
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

#ifndef __TFBAYES_DPM_DPM_TFBS_OPTIONS_HH__
#define __TFBAYES_DPM_DPM_TFBS_OPTIONS_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <string>
#include <vector>
#include <ostream>

#include <tfbayes/dpm/dpm-partition.hh>

typedef std::matrix<double>             baseline_lengths_t;
typedef std::list<std::string>          baseline_names_t;
typedef std::list<std::matrix<double> > baseline_priors_t;
typedef std::vector<double>             baseline_weights_t;

typedef struct {
        std::string phylogenetic_file;
        std::string alignment_file;
        double alpha;
        double discount;
        double lambda;
        double initial_temperature;
        bool   block_samples;
        size_t block_samples_period;
        size_t metropolis_proposals;
        bool   optimize;
        size_t optimize_period;
        std::string process_prior;
        std::string background_model;
        std::matrix<double> background_alpha;
        size_t background_context;
        std::vector<double> background_beta;
        std::vector<double> background_gamma;
        std::string background_cache;
        std::vector<double> background_weights;
        baseline_lengths_t baseline_lengths;
        baseline_names_t   baseline_names;
        baseline_priors_t  baseline_priors;
        baseline_weights_t baseline_weights;
        size_t population_size;
        size_t threads;
        std::string socket_file;
        size_t verbose;
} tfbs_options_t;

std::ostream& operator<<(std::ostream& o, const tfbs_options_t& options);

#endif /* __TFBAYES_DPM_DPM_TFBS_OPTIONS_HH__ */
