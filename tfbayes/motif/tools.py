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

import math

from ..uipac.alphabet import DNA

# config tools
# ------------------------------------------------------------------------------

def read_vector(config, section, option, converter):
    vector_str = config.get(section, option)
    vector     = map(converter, vector_str.split(' '))
    return vector

def read_matrix(config, section, option, converter):
    matrix_str = config.get(section, option)
    matrix     = []
    for line in matrix_str.split('\n'):
        if line != '':
            matrix.append([converter(a) for a in line.split(' ')])
    return matrix

# counts -> frequencies
# ------------------------------------------------------------------------------

def compute_frequencies(counts):
    sums = [ sum(map(lambda m: m[j], counts))     for j in range(0, len(counts[0])) ]
    return [ [ float(counts[i][j])/float(sums[j]) for j in range(0, len(counts[0])) ] for i in range(len(counts)) ]

# frequencies -> pwm
# ------------------------------------------------------------------------------

def compute_pwm(config_parser, cluster_name):
    if config_parser.has_option('Cluster', 'cluster_bg'):
        bg_freq   = read_matrix(config_parser, 'Cluster', 'cluster_bg',  float)
    else:
        bg_freq = [[0.0], [0.0], [0.0], [0.0]]
        bg_freq[DNA.code('A')][0] = 0.3
        bg_freq[DNA.code('C')][0] = 0.2
        bg_freq[DNA.code('G')][0] = 0.2
        bg_freq[DNA.code('T')][0] = 0.3
    tfbs_counts = read_matrix(config_parser, 'Cluster', cluster_name, float)
    tfbs_counts = [ [ tfbs_counts[i][j]+1 for j in range(len(tfbs_counts[0])) ] for i in range(4) ]
    tfbs_freq   = compute_frequencies(tfbs_counts)
    return [ [ math.log(tfbs_freq[i][j]/bg_freq[i][0], 2) for j in range(len(tfbs_freq[0])) ] for i in range(4) ]

# score of a sequence given a pwm
# ------------------------------------------------------------------------------

def compute_score(pwm, sequence):
    result = 0
    for j in range(len(sequence)):
        if sequence[j] == 'A' or sequence[j] == 'a':
            code = DNA.code('A')
            result += pwm[code][j]
        if sequence[j] == 'C' or sequence[j] == 'c':
            code = DNA.code('C')
            result += pwm[code][j]
        if sequence[j] == 'G' or sequence[j] == 'g':
            code = DNA.code('G')
            result += pwm[code][j]
        if sequence[j] == 'T' or sequence[j] == 't':
            code = DNA.code('T')
            result += pwm[code][j]
    return result

# reverse complement
# ------------------------------------------------------------------------------

def reverse_complement(motif):
    revcomp = [ [ 0.0 for i in range(len(motif[0])) ] for j in range(len(motif)) ]
    for i in range(len(motif[0])):
        revcomp[DNA.code('T')][-i-1] = motif[DNA.code('A')][i]
        revcomp[DNA.code('G')][-i-1] = motif[DNA.code('C')][i]
        revcomp[DNA.code('C')][-i-1] = motif[DNA.code('G')][i]
        revcomp[DNA.code('A')][-i-1] = motif[DNA.code('T')][i]
    return revcomp
