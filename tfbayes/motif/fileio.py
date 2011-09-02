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

# config parser utilities
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

def write_vector(config, section, option, vector):
    config.set(section, option, " ".join(map(str, vector)))

def write_matrix(config, section, option, matrix):
    config.set(section, option, "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), matrix)))

# load and save motifs
# ------------------------------------------------------------------------------

def save_motifs(config_parser, bg, motifs, components):
    if not config_parser.has_section('Cluster'):
        config_parser.add_section('Cluster')
    cluster_names = []
    for n, motif, comp in zip(range(0, len(motifs)), motifs, components):
        cluster_name = 'cluster_%d' % n
        cluster_names.append(cluster_name)
        write_matrix(config_parser, 'Cluster', cluster_name, motif)
        config_parser.set('Cluster', '%s_components' % cluster_name, str(comp))
    write_vector(config_parser, 'Cluster', 'cluster', cluster_names)
    write_matrix(config_parser, 'Cluster', 'cluster_bg', bg)

def load_motifs(config_parser):
    bg = read_matrix(config_parser, 'Cluster', 'cluster_bg', float)
    motifs = []
    components = []
    cluster_names = read_vector(config_parser, 'Cluster', 'cluster', str)
    for cluster_name in cluster_names:
        motifs.append(read_matrix(config_parser, 'Cluster', cluster_name, float))
        components.append(int(config_parser.get('Cluster', '%s_components' % cluster_name)))
    return bg, motifs, components
