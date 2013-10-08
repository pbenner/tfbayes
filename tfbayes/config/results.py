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

import ConfigParser

from tools import *

from tfbayes.cluster   import parse_partition_list
from tfbayes.interface import *
from tfbayes.dpm.dpm_tfbs_interface import *

# default results config
# ------------------------------------------------------------------------------

def default_results_config():
    results_config = sampling_history_t()
    return results_config

# parse results config
# ------------------------------------------------------------------------------

def parse_results_config(config_file, results_config):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)
    if not config_parser.has_section('Result'):
        raise IOError("Invalid configuration file.")
    if config_parser.has_option('Result', 'components'):
        results_config.components = read_matrix(config_parser, 'Result', 'components', int)
    if config_parser.has_option('Result', 'likelihood'):
        results_config.likelihood = read_matrix(config_parser, 'Result', 'likelihood', float)
    if config_parser.has_option('Result', 'posterior'):
        results_config.posterior = read_matrix(config_parser, 'Result', 'posterior', float)
    if config_parser.has_option('Result', 'temperature'):
        results_config.temperature = read_matrix(config_parser, 'Result', 'temperature', float)
    if config_parser.has_option('Result', 'switches'):
        results_config.switches = read_matrix(config_parser, 'Result', 'switches', float)
    if config_parser.has_option('Result', 'map_partition'):
        results_config.map_partition = config_parser.get('Result', 'map_partition')
    if config_parser.has_option('Result', 'mean_partition'):
        results_config.mean_partition = config_parser.get('Result', 'mean_partition')
    if config_parser.has_option('Result', 'median_partition'):
        results_config.median_partition = config_parser.get('Result', 'median_partition')
    if config_parser.has_option('Result', 'partitions'):
        results_config.partitions = parse_partition_list(config_parser.get('Result', 'partitions'))

# save results config
# ------------------------------------------------------------------------------

def save_results_config(config_file, results_config):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.add_section('Result')
    write_matrix(config_parser, 'Result', 'components',  results_config.components)
    write_matrix(config_parser, 'Result', 'likelihood',  results_config.likelihood)
    write_matrix(config_parser, 'Result', 'posterior',   results_config.posterior)
    write_matrix(config_parser, 'Result', 'switches',    results_config.switches)
    write_matrix(config_parser, 'Result', 'temperature', results_config.temperature)
    if results_config.has_key('partitions'):
        config_parser.set('Result', 'partitions',       results_config.partitions)
    if results_config.has_key('map_partition'):
        config_parser.set('Result', 'map_partition',    results_config.map_partition)
    if results_config.has_key('mean_partition'):
        config_parser.set('Result', 'mean_partition',   results_config.mean_partition)
    if results_config.has_key('median_partition'):
        config_parser.set('Result', 'median_partition', results_config.median_partition)
    with open(config_file, 'wb') as config_fp:
        config_parser.write(config_fp)
