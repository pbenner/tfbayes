/* Copyright (C) 2015 Philipp Benner
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

#ifndef __TFBAYES_CONFIG_PARTITION_MAIN_HH__
#define __TFBAYES_CONFIG_PARTITION_MAIN_HH__

#include <string>

#include <tfbayes/dpm/dpm-partition.hh>

dpm_partition_list_t parse_partition_list(FILE * file);
dpm_partition_list_t parse_partition_list(const std::string& filename);
dpm_partition_list_t parse_partition_list_str(const std::string& str);
dpm_partition_t parse_partition_str(const std::string& str);

#endif /* __TFBAYES_CONFIG_PARTITION_MAIN_HH__ */
