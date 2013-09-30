#! /usr/bin/env python

# Copyright (C) 2011, 2012 Philipp Benner
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

from cluster import cluster_t

# parse results partition
# ------------------------------------------------------------------------------

class dpm_subset_t():
    # the identifier gives the cluster type, i.e. which component model
    # or prior the subset represents
    identifier = None
    subset     = None
    def __str__(self):
        result = "%s:{" % self.identifier
        for i, index in enumerate(self.subset):
            result += str(index)
            if not i is len(self.subset)-1:
                result += ", "
        result += "}"
        return result

class dpm_partition_t(list):
    def __str__(self):
        result = ""
        for i, dpm_subset in enumerate(self):
            result += str(dpm_subset)
            if not i is len(self)-1:
                result += ", "
        return result

def parse_partition_identifier(partition_str, end):
    # start is the position of '{', we want the identifier
    # before this bracket
    start = 0
    if not partition_str[end-1] == ':':
        raise ValueError('Invalid MAP partition.')
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
            start = stack.pop()
            element = subset_str[start + 1: i].split(',')
            yield (int(element[0]), int(element[1]))

def parse_partition(partition_str):
    partition = dpm_partition_t()
    for identifier, subset in parse_partition_subsets(partition_str):
        dpm_subset = dpm_subset_t()
        dpm_subset.identifier = identifier
        dpm_subset.subset     = list(parse_partition_elements(subset))
        partition.append(dpm_subset)
    return partition

# parse partition list
# ------------------------------------------------------------------------------

def parse_partition_list(partition_list_str):
    partition_list = []
    for partition_str in partition_list_str.split('\n'):
        partition = parse_partition(partition_str)
        # there might be empty lines, filter them out...
        if partition:
            partition_list.append(partition)
    return partition_list

# sort cluster
# ------------------------------------------------------------------------------

def update_identifier(cluster_list):
    for idx, cluster in enumerate(cluster_list):
        cluster.identifier = idx
    return cluster_list

def permutation_indices(scores):
     return sorted(range(len(scores)), key = scores.__getitem__, reverse=True)

def sort_cluster_list(cluster_list):
    # the score is the criterion according to which the list of motifs
    # is sorted
    scores = [ cluster.average_counts() for cluster in cluster_list ]
    # sort this list and obtain the permutation pi
    pi     = permutation_indices(scores)
    # return a sorted list of motifs
    return update_identifier([ cluster_list[i] for i in pi ])

# generate motifs
# ------------------------------------------------------------------------------

def generate_cluster(sequences, dpm_subset, identifier, sampler_config):
    counts       = [ [ 0.0 for j in range(sampler_config['tfbs_length']) ] for i in range(4) ]
    counts_gap   = [ 0.0 ] * sampler_config['tfbs_length']
    alpha        = sampler_config['baseline_priors'][dpm_subset.identifier][0:4]
    alpha_gap    = sampler_config['baseline_priors'][dpm_subset.identifier][4]
    components   = len(dpm_subset.subset)
    cluster_type = dpm_subset.identifier
    for position in dpm_subset.subset:
        s = position[0]
        p = position[1]
        # loop over the motif
        for j in range(sampler_config['tfbs_length']):
            # loop over all nucleotides plus counts for gaps
            for i in range(4):
                counts[i][j] += sequences[s][p+j][i]
            counts_gap[j] += sequences[s][p+j][4]
    return cluster_t(counts, counts_gap, alpha, alpha_gap, components, identifier, cluster_type)

def generate_cluster_list(sequences, sampler_config, results_config, which = 'map'):
    """Loop through the partition and for each subset generate a motif."""
    cluster_list = []
    # receive partition from results file
    if which == 'map':
        if not results_config.has_key('map_partition'):
            raise IOError("MAP partition is not available in results file.")
        partition = parse_partition(results_config['map_partition'])
    elif which == 'mean':
        if not results_config.has_key('mean_partition'):
            raise IOError("Mean partition is not available in results file.")
        partition = parse_partition(results_config['mean_partition'])
    elif which == 'median':
        if not results_config.has_key('median_partition'):
            raise IOError("Median partition is not available in results file.")
        partition = parse_partition(results_config['median_partition'])
    else:
        raise IOError("`%s' is not a valid partition type." % which)
    # generate cluster list
    for idx, dpm_subset in enumerate(partition):
        cluster_list.append(generate_cluster(sequences, dpm_subset, idx, sampler_config))
    return sort_cluster_list(cluster_list)
