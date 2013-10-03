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

#ifndef DATATYPES_HH
#define DATATYPES_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <vector>
#include <string>

#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/linalg.hh>
#include <tfbayes/dpm/index.hh>

// cluster structures
////////////////////////////////////////////////////////////////////////////////

typedef ssize_t cluster_tag_t;
typedef std::string baseline_tag_t;

typedef enum {
        cluster_event_empty, cluster_event_nonempty,
        cluster_event_add_word, cluster_event_remove_word
} cluster_event_t;

#include <tfbayes/dpm/dpm-partition.hh>

class sampling_history_t {
public:
        std::matrix<double> switches;
        std::matrix<double> likelihood;
        std::matrix<double> posterior;
        std::matrix<double> components;
        std::matrix<double> temperature;
        std::vector<dpm_partition_t> partitions;
};

#endif /* DATATYPES_HH */
