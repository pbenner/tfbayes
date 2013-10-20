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

from ..uipac import *
from ..config.tools import read_vector, write_vector

# pwm to string
# ------------------------------------------------------------------------------

def vector_to_string(vector):
    string = ''
    for scalar in vector:
        string += '%6.2f ' % scalar
    return string

def pwm_to_string(pwm):
    string = ''
    for idx, vector in enumerate(pwm):
        string += '%s: ' % DNA.decode(idx),
        string += vector_to_string(vector)
        string += '\n'
    return string

# counts to string
# ------------------------------------------------------------------------------

def counts_to_string(counts, counts_gaps):
    string = ''
    for idx, vector in enumerate(counts+[counts_gaps]):
        string += '%s: ' % DNA.decode(idx)
        string += vector_to_string(vector)
        string += '\n'
    return string

# print pwm
# ------------------------------------------------------------------------------

def print_pwm(pwm):
    print pwm_to_string(pwm)

# print counts
# ------------------------------------------------------------------------------

def print_counts(counts, counts_gaps):
    print counts_to_string(counts, counts_gaps)

# parse pwm
# ------------------------------------------------------------------------------

def parse_pwm_line(line):
    m = re.match('([ACGTacgt]):', line[0])
    if not m:
        raise IOError("Parsing of PWM failed!")
    nucleotide = m.group(1)
    entry      = map(float, line[1:])
    return nucleotide, entry

def parse_pwm(pwm_string):
    pwm    = [[]] * 4
    for line in pwm_string.strip().split('\n'):
        line = filter(lambda x: not x is '', line.strip().split(' '))
        if not len(line) is 0:
            nucleotide, entry = parse_pwm_line(line)
            pwm[DNA.code(nucleotide)] = entry
    return pwm

# parse counts
# ------------------------------------------------------------------------------

def parse_counts_line(line):
    m = re.match('([ACGTacgt-]):', line[0])
    if not m:
        raise IOError("Parsing of COUNTS failed!")
    nucleotide = m.group(1)
    entry      = map(float, line[1:])
    return nucleotide, entry

def parse_counts(counts_string):
    counts    = [[]] * 5
    for line in counts_string.strip().split('\n'):
        line = filter(lambda x: not x is '', line.strip().split(' '))
        if not len(line) is 0:
            nucleotide, entry = parse_counts_line(line)
            counts[DNA.code(nucleotide)] = entry
    return counts[:4], counts[4]

# read counts from config file
# ------------------------------------------------------------------------------

def read_counts(cluster_parser, cluster_name):
    counts_string = cluster_parser.get('Cluster', cluster_name)
    return parse_counts(counts_string)

# load and save cluster
# ------------------------------------------------------------------------------

def load_cluster(cluster_parser, sampler_config, cluster_name):
    result       = re.match('cluster_([0-9]+)', cluster_name)
    counts, counts_gap = read_counts(cluster_parser, cluster_name)
    components   = int(cluster_parser.get('Cluster', '%s_components' % cluster_name))
    identifier   = int(result.group(1))
    cluster_type = cluster_parser.get('Cluster', '%s_type' % cluster_name)
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
    # write cluster counts and alpha
    cluster_parser.set('Cluster', cluster_name, '\n'+counts_to_string(cluster.counts, cluster.counts_gap).strip())
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
