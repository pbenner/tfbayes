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

import sys

from cluster import cluster_t
from ..uipac import DNA

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

def get_baseline_prior(id, sampler_config):
    for baseline_prior, baseline_name in zip(sampler_config.baseline_priors,
                                             sampler_config.baseline_names):
        if id.name == baseline_name:
            return map(list, zip(*baseline_prior))
    raise IOError("Baseline prior not found!")

def generate_cluster(sequences, dpm_subset, idx, sampler_config, index_error = True):
    length         = dpm_subset.model_id().length
    counts         = [ [ 0.0 for j in range(length) ] for i in range(4) ]
    counts_gap     = [ 0.0 ] * length
    baseline_prior = get_baseline_prior(dpm_subset.model_id(), sampler_config)
    alpha          = [ [ baseline_prior[i][0] for j in range(length) ] for i in range(4) ]
    alpha_gap      = baseline_prior[4] * length
    components     = len(dpm_subset)
    for r in dpm_subset:
        s = r.index()[0]
        p = r.index()[1]
        # loop over the motif
        for j in range(r.length()):
            # loop over all nucleotides plus counts for gaps
            if r.reverse():
                try: sequences[s][p+j][DNA.complement(i)]
                except IndexError:
                    if index_error:
                        raise
                    else:
                        print >> sys.stderr, "Warning: index %d:%d out of range!" % (s, p+j)
                        continue
                for i in range(4):
                    counts[i][r.length()-j-1] += sequences[s][p+j][DNA.complement(i)]
            else:
                try: sequences[s][p+j][i]
                except IndexError:
                    if index_error:
                        raise
                    else:
                        print >> sys.stderr, "Warning: index %d:%d out of range!" % (s, p+j)
                        continue
                for i in range(4):
                    counts[i][j] += sequences[s][p+j][i]
            counts_gap[j] += sequences[s][p+j][4]
    return cluster_t(counts, counts_gap, alpha, alpha_gap, components, idx, dpm_subset.model_id(), sites=dpm_subset)

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
