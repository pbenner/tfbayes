#! /usr/bin/env python

# Copyright (C) 2011 Philipp Benner
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

from ..config import *

# load and save cluster
# ------------------------------------------------------------------------------

def load_cluster(config_parser, cluster_name):
    result      = re.match('cluster_([0-9]+)', cluster_name)
    counts      = read_matrix(config_parser, 'Cluster', cluster_name, float)[:4]
    counts_gap  = read_matrix(config_parser, 'Cluster', cluster_name, float)[4]
    components  = int(config_parser.get('Cluster', '%s_components' % cluster_name))
    identifier  = int(result.group(1))
    return cluster_t(counts, counts_gap, components, identifier)

def load_cluster_list(config_parser):
    cluster_list  = []
    cluster_names = read_vector(config_parser, 'Cluster', 'cluster', str)
    for cluster_name in cluster_names:
        cluster_list.append(load_cluster(config_parser, cluster_name))
    return cluster_list

def save_cluster(config_parser, cluster):
    if not config_parser.has_section('Cluster'):
        config_parser.add_section('Cluster')
    # name of the cluster in the config file
    cluster_name = 'cluster_%d' % cluster.identifier
    # write cluster counts
    write_matrix(config_parser, 'Cluster', cluster_name, cluster.counts)
    # write number of cluster components
    config_parser.set('Cluster', '%s_components' % cluster_name, cluster.components)
    return cluster_name

def save_cluster_list(config_parser, cluster_list):
    cluster_names = []
    for cluster in cluster_list:
        cluster_name = save_cluster(config_parser, cluster)
        cluster_names.append(cluster_name)
    write_vector(config_parser, 'Cluster', 'cluster', cluster_names)
