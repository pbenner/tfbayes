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

typedef std::list<std::matrix<double> > baseline_priors_t;
typedef std::list<std::string> baseline_tags_t;

typedef struct {
        std::string phylogenetic_file;
        std::string alignment_file;
        size_t tfbs_length;
        double alpha;
        double discount;
        double lambda;
        double initial_temperature;
        std::string process_prior;
        std::string background_model;
        std::matrix<double> background_alpha;
        size_t background_context;
        std::vector<double> background_gamma;
        std::string background_cache;
        std::string background_weights;
        std::vector<double> baseline_weights;
        baseline_priors_t baseline_priors;
        baseline_tags_t baseline_tags;
        size_t population_size;
        size_t threads;
        std::string socket_file;
} tfbs_options_t;

std::ostream& operator<<(std::ostream& o, const tfbs_options_t& options);

#endif /* __TFBAYES_DPM_DPM_TFBS_OPTIONS_HH__ */
