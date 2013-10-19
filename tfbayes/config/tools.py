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

def read_vector(config_parser, section, option, converter):
    vector_str = config_parser.get(section, option).replace('  ',' ').strip()
    vector     = map(converter, vector_str.split(' '))
    return vector

def read_matrix(config_parser, section, option, converter):
    matrix_str = config_parser.get(section, option)
    matrix     = []
    for line in matrix_str.split('\n'):
        line = line.replace('  ',' ').strip()
        if line != '':
            matrix.append([converter(a) for a in line.split(' ')])
    return matrix

def write_vector(config, section, option, vector):
    config.set(section, option, " ".join(map(str, vector)))

def write_matrix(config, section, option, matrix, converter = None):
    if converter:
        config.set(section, option, "\n"+"\n".join(map(lambda arg: " ".join(map(lambda a: str(int(a)), arg)), matrix)))
    else:
        config.set(section, option, "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), matrix)))
