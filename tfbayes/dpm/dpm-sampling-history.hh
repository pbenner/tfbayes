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

#ifndef __TFBAYES_DPM_SAMPLING_HISTORY_HH__
#define __TFBAYES_DPM_SAMPLING_HISTORY_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/utility/linalg.hh>
#include <tfbayes/dpm/dpm-partition.hh>

class sampling_history_t {
public:
        std::matrix<double>  switches;
        std::matrix<double>  likelihood;
        std::matrix<double>  posterior;
        std::matrix<double>  components;
        std::matrix<double>  temperature;
        std::matrix<double>  cluster_sizes;
        dpm_partition_list_t partitions;
};

#endif /* __TFBAYES_DPM_SAMPLING_HISTORY_HH__ */
