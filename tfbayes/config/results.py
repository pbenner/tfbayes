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

import ConfigParser

from tools import *

# parse results config
# ------------------------------------------------------------------------------

def parse_results_config(config_file, results_config):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)
    if not config_parser.has_section('Result'):
        raise IOError("Invalid configuration file.")
    if config_parser.has_option('Result', 'components'):
        results_config['components'] = read_matrix(config_parser, 'Result', 'components', int)
    if config_parser.has_option('Result', 'likelihood'):
        results_config['likelihood'] = read_matrix(config_parser, 'Result', 'likelihood', float)
    if config_parser.has_option('Result', 'posterior'):
        results_config['posterior'] = read_matrix(config_parser, 'Result', 'posterior', float)
    if config_parser.has_option('Result', 'switches'):
        results_config['switches'] = read_matrix(config_parser, 'Result', 'switches', float)
    if config_parser.has_option('Result', 'graph'):
        results_config['graph'] = config_parser.get('Result', 'graph')
    if config_parser.has_option('Result', 'map_partition'):
        results_config['map_partition'] = config_parser.get('Result', 'map_partition')

# save results config
# ------------------------------------------------------------------------------

def save_results_config(config_file, results_config):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.add_section('Result')
    write_matrix(config_parser, 'Result', 'components', results_config['components'])
    write_matrix(config_parser, 'Result', 'likelihood', results_config['likelihood'])
    write_matrix(config_parser, 'Result', 'posterior',  results_config['posterior'])
    write_matrix(config_parser, 'Result', 'switches',   results_config['switches'])
    config_parser.set('Result', 'graph', results_config['graph'])
    config_parser.set('Result', 'map_partition', results_config['map_partition'])
    with open(config_file, 'wb') as config_fp:
        config_parser.write(config_fp)
