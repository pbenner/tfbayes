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

#ifndef DPM_TFBS_OPTIONS_HH
#define DPM_TFBS_OPTIONS_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <vector>

#include <tfbayes/dpm/dpm-partition.hh>

typedef struct {
        size_t tfbs_length;
        double alpha;
        double discount;
        double lambda;
        double initial_temperature;
        bool construct_graph;
        std::string process_prior;
        std::string background_model;
        std::matrix<double> background_alpha;
        size_t background_context;
        std::string background_weights;
        std::vector<double> baseline_weights;
        std::vector<std::matrix<double> > baseline_priors;
        std::vector<std::string> baseline_tags;
        dpm_partition_t partition;
        std::string socket_file;
} tfbs_options_t;

#endif /* DPM_TFBS_OPTIONS_HH */
