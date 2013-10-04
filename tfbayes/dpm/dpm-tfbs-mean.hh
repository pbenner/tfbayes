/* Copyright (C) 2013 Philipp Benner
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

#ifndef DPM_TFBS_MEDIAN_HH
#define DPM_TFBS_MEDIAN_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <tfbayes/dpm/datatypes.hh>
#include <tfbayes/dpm/dpm-partition.hh>

size_t
dpm_tfbs_mean(const std::vector<dpm_partition_t>& partitions,
              const std::vector<size_t>& sizes, size_t tfbs_length,
              bool verbose = false);

size_t
dpm_tfbs_median(const std::vector<dpm_partition_t>& partitions,
                const std::vector<size_t>& sizes, size_t tfbs_length,
                bool verbose = false);

#endif /* DPM_TFBS_MEDIAN_HH */
