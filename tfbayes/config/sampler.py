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

# parse config
# ------------------------------------------------------------------------------

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def parse_sampler_config(config_file, options):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)
    if not config_parser.has_section('TFBS-Sampler'):
        raise IOError("Invalid configuration file.")
    options['seq_file'] = config_parser.get('TFBS-Sampler', 'sequences')
    if config_parser.has_option('TFBS-Sampler', 'alpha'):
        options['alpha'] = float(config_parser.get('TFBS-Sampler', 'alpha'))
    if config_parser.has_option('TFBS-Sampler', 'discount'):
        options['discount'] = float(config_parser.get('TFBS-Sampler', 'discount'))
    if config_parser.has_option('TFBS-Sampler', 'lambda'):
        options['lambda'] = float(config_parser.get('TFBS-Sampler', 'lambda'))
    if config_parser.has_option('TFBS-Sampler', 'process-prior'):
        options['process_prior'] = config_parser.get('TFBS-Sampler', 'process-prior').strip()
    if config_parser.has_option('TFBS-Sampler', 'background-model'):
        options['background_model'] = config_parser.get('TFBS-Sampler', 'background-model').strip()
    if config_parser.has_option('TFBS-Sampler', 'background-alpha'):
        options['background_alpha'] = read_matrix(config_parser, 'TFBS-Sampler', 'background-alpha', float)
    if config_parser.has_option('TFBS-Sampler', 'background-context'):
        options['background_context'] = int(config_parser.get('TFBS-Sampler', 'background-context'))
    if config_parser.has_option('TFBS-Sampler', 'background-weights'):
        options['background_weights'] = config_parser.get('TFBS-Sampler', 'background-weights')
    if config_parser.has_option('TFBS-Sampler', 'tfbs-length'):
        options['tfbs_length'] = int(config_parser.get('TFBS-Sampler', 'tfbs-length'))
    if config_parser.has_option('TFBS-Sampler', 'construct-graph'):
        options['construct_graph'] = str2bool(config_parser.get('TFBS-Sampler', 'construct_graph'))
    if config_parser.has_option('TFBS-Sampler', 'metropolis-optimize'):
        options['metropolis_optimize'] = str2bool(config_parser.get('TFBS-Sampler', 'metropolis-optimize'))
    if config_parser.has_option('TFBS-Sampler', 'population-size'):
        options['population_size'] = int(config_parser.get('TFBS-Sampler', 'population-size'))
    if config_parser.has_option('TFBS-Sampler', 'save'):
        options['save'] = config_parser.get('TFBS-Sampler', 'save')
    if config_parser.has_option('TFBS-Sampler', 'samples'):
        options['samples'] = tuple(map(int, config_parser.get('TFBS-Sampler', 'samples').split(":")))
    if config_parser.has_option('TFBS-Sampler', 'baseline-weights'):
        options['baseline_weights'] = read_vector(config_parser, 'TFBS-Sampler', 'baseline-weights', int)
    if config_parser.has_option('TFBS-Sampler', 'baseline-priors'):
        options['baseline_priors'] = []
        for prior_name in read_vector(config_parser, 'TFBS-Sampler', 'baseline-priors', str):
            prior = read_matrix(config_parser, 'TFBS-Sampler', prior_name, float)
            options['baseline_priors'].append(prior)
    if config_parser.has_option('TFBS-Sampler', 'socket-file'):
        options['socket_file'] = config_parser.get('TFBS-Sampler', 'socket-file').strip()
