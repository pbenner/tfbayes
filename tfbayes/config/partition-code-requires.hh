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

#ifndef __TFBAYES_CONFIG_PARTITION_CODE_HH__
#define __TFBAYES_CONFIG_PARTITION_CODE_HH__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <new>
#include <stack>
#include <string>
#include <cstdlib>

#include <tfbayes/dpm/dpm-partition.hh>

struct node_t {
	dpm_partition_list_t partition_list;
	dpm_partition_t partition;
        dpm_subset_t subset;
        range_t range;
        std::string string;
        size_t integer;
};

#define YYSTYPE std::stack<node_t> *

typedef void* yyscan_t;
struct context_t {
        dpm_partition_list_t partition_list;
	yyscan_t scanner;
};

#endif /* __TFBAYES_CONFIG_PARTITION_CODE_HH__ */
