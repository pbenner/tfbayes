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

from ..uipac  import *
from ..config import *
from cluster  import *

# load and save cluster
# ------------------------------------------------------------------------------

def load_cluster(cluster_parser, sampler_config, cluster_name):
    result       = re.match('cluster_([0-9]+)', cluster_name)
    counts       = read_matrix(cluster_parser, 'Cluster', cluster_name, float)[:4]
    counts_gap   = read_matrix(cluster_parser, 'Cluster', cluster_name, float)[4]
    components   = int(cluster_parser.get('Cluster', '%s_components' % cluster_name))
    identifier   = int(result.group(1))
    cluster_type = int(cluster_parser.get('Cluster', '%s_type' % cluster_name))
    alpha        = sampler_config['baseline_priors'][cluster_type][0:4]
    alpha_gap    = sampler_config['baseline_priors'][cluster_type][4]
    return cluster_t(counts, counts_gap, alpha, alpha_gap, components, identifier, cluster_type)

def load_cluster_list(cluster_parser, sampler_config):
    cluster_list  = []
    cluster_names = read_vector(cluster_parser, 'Cluster', 'cluster', str)
    for cluster_name in cluster_names:
        cluster_list.append(load_cluster(cluster_parser, sampler_config, cluster_name))
    return cluster_list

def save_cluster(cluster_parser, cluster):
    if not cluster_parser.has_section('Cluster'):
        cluster_parser.add_section('Cluster')
    # name of the cluster in the config file
    cluster_name = 'cluster_%d' % cluster.identifier
    # write cluster counts
    write_matrix(cluster_parser, 'Cluster', cluster_name, cluster.counts+[cluster.counts_gap])
    # write number of cluster components
    cluster_parser.set('Cluster', '%s_components' % cluster_name, cluster.components)
    # write cluster type (i.e. which prior was used)
    cluster_parser.set('Cluster', '%s_type' % cluster_name, cluster.cluster_type)
    return cluster_name

def save_cluster_list(cluster_parser, cluster_list):
    cluster_names = []
    for cluster in cluster_list:
        cluster_name = save_cluster(cluster_parser, cluster)
        cluster_names.append(cluster_name)
    write_vector(cluster_parser, 'Cluster', 'cluster', cluster_names)

# print pwm
# ------------------------------------------------------------------------------

def print_vector(vector):
    for scalar in vector:
        print '%5.2f ' % scalar,
    print

def print_pwm(matrix):
    for idx, vector in enumerate(matrix):
        print '%s: ' % DNA.decode(idx),
        print_vector(vector)

# parse pwm
# ------------------------------------------------------------------------------

def parse_pwm_line(line):
    m = re.match('([ACGTacgt]):', line[0])
    if not m:
        raise ValueError("Parsing of PWM failed!")
    nucleotide = m.group(1)
    entry      = map(float, line[1:])
    return nucleotide, entry

def parse_pwm(pwm_file):
    pwm_fp = open(pwm_file, 'r')
    pwm    = [[]] * 4
    for line in pwm_fp:
        line = filter(lambda x: not x is '', line.strip().split(' '))
        if not len(line) is 0:
            nucleotide, entry = parse_pwm_line(line)
            pwm[DNA.code(nucleotide)] = entry
    pwm_fp.close()
    return pwm
