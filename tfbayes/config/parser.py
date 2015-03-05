#! /usr/bin/env python

# Copyright (C) 2011-2013 Philipp Benner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import re

from ..interface import *
from ..dpm       import *

# parse results partition
# ------------------------------------------------------------------------------

def parse_partition_identifier(partition_str, end):
    # start is the position of '{', we want the identifier
    # before this bracket
    start = 0
    if not partition_str[end-1] == ':':
        raise IOError('Invalid partition.')
    for i in range(end-1,-1,-1):
        if (partition_str[i] == ' ' or partition_str[i] == ','):
            start = i+1
            break
    return partition_str[start:end-1]

def parse_partition_subsets(partition_str):
    """Split the partition string into a list of its subsets."""
    stack = []
    for i, c in enumerate(partition_str):
        if c == '{':
            stack.append(i)
        elif c == '}' and stack:
            start      = stack.pop()
            identifier = parse_partition_identifier(partition_str, start)
            yield identifier, partition_str[start + 1: i]

def parse_partition_elements(subset_str):
    """Split the string of a subset into its elements and parse them."""
    stack = []
    for i, c in enumerate(subset_str):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            start  = stack.pop()
            index  = subset_str[start + 1: i].split(',')
            m      = re.match(":([0-9]+)(!?)", subset_str[i+1:len(subset_str)])
            if m == None:
                raise IOError('Invalid partition.')
            length = int(m.groups()[0])
            rc     = m.groups()[1] == "!"
            yield range_t(seq_index_t(int(index[0]), int(index[1])), length, rc)

def parse_partition(partition_str):
    partition = dpm_partition_t()
    for identifier, subset in parse_partition_subsets(partition_str):
        model_id        = model_id_t()
        model_id.name   = identifier.split(":")[0]
        model_id.length = int(identifier.split(":")[1])
        dpm_subset = dpm_subset_t(model_id)
        for r in parse_partition_elements(subset):
            dpm_subset.insert(r)
        partition.append(dpm_subset)
    return partition

# parse partition list
# ------------------------------------------------------------------------------

def parse_partition_list(partition_list_str):
    partition_list = dpm_partition_list_t()
    for partition_str in partition_list_str.split('\n'):
        partition = parse_partition(partition_str)
        # there might be empty lines, filter them out...
        if partition:
            partition_list.append(partition)
    return partition_list
