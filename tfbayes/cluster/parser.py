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

from cluster import cluster_t

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

def get_baseline_prior(tag, sampler_config):
    for baseline_prior, baseline_tag in zip(sampler_config.baseline_priors,
                                            sampler_config.baseline_tags):
        if tag == baseline_tag:
            return map(list, zip(*baseline_prior))
    raise IOError("Baseline prior not found!")

def generate_cluster(sequences, dpm_subset, idx, sampler_config):
    counts         = [ [ 0.0 for j in range(sampler_config.tfbs_length) ] for i in range(4) ]
    counts_gap     = [ 0.0 ] * sampler_config.tfbs_length
    baseline_prior = get_baseline_prior(dpm_subset.dpm_subset_tag(), sampler_config)
    alpha          = baseline_prior[0:4]
    alpha_gap      = baseline_prior[4]
    components     = len(dpm_subset)
    for position in dpm_subset:
        s = position[0]
        p = position[1]
        # loop over the motif
        for j in range(sampler_config.tfbs_length):
            # loop over all nucleotides plus counts for gaps
            for i in range(4):
                counts[i][j] += sequences[s][p+j][i]
            counts_gap[j] += sequences[s][p+j][4]
    return cluster_t(counts, counts_gap, alpha, alpha_gap, components, idx, dpm_subset.dpm_subset_tag())

def generate_cluster_list(sequences, sampler_config, results_config, which = 'map'):
    """Loop through the partition and for each subset generate a motif."""
    cluster_list = []
    # receive partition from results file
    if which == 'map':
        if not results_config.has_key('map_partition') or not results_config['map_partition']:
            raise IOError("MAP partition is not available in results file. Use `tfbayes-estimate' to compute it.")
        partition = results_config['map_partition']
    elif which == 'mean':
        if not results_config.has_key('mean_partition') or not results_config['mean_partition']:
            raise IOError("Mean partition is not available in results file. Use `tfbayes-estimate' to compute it.")
        partition = results_config['mean_partition']
    elif which == 'median':
        if not results_config.has_key('median_partition') or not results_config['median_partition']:
            raise IOError("Median partition is not available in results file. Use `tfbayes-estimate' to compute it.")
        partition = results_config['median_partition']
    else:
        raise IOError("`%s' is not a valid partition type." % which)
    # generate cluster list
    for idx, dpm_subset in enumerate(partition):
        cluster_list.append(generate_cluster(sequences, dpm_subset, idx, sampler_config))
    return sort_cluster_list(cluster_list)
