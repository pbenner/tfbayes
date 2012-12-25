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

import copy
import ConfigParser

from tools import *

# default sampler config
# ------------------------------------------------------------------------------

default_sampler_config_ = {
    'alpha'               : 0.05,
    'discount'            : 0.0,
    'lambda'              : 0.01,
    'construct_graph'     : False,
    'process_prior'       : "pitman-yor process",
    'background_model'    : "independence-dirichlet",
    'background_alpha'    : [[1],[1],[1],[1],[5]],
    'background_context'  : 2,
    'background_weights'  : 'decay',
    'population_size'     : 1,
    'tfbs_length'         : 10,
    'seq_file'            : None,
    'baseline_names'      : [],
    'baseline_weights'    : {},
    'baseline_priors'     : {},
    'socket_file'         : '',
    'samples'             : (1000,100),
    'save'                : None
    }

def default_sampler_config():
    return copy.deepcopy(default_sampler_config_)

# parse config
# ------------------------------------------------------------------------------

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def normalize_dict(d):
    norm = float(sum(d.viewvalues()))
    for k, v in d.iteritems():
        d[k] = d[k] / norm

def generate_baseline(sampler_config):
    length = sampler_config['tfbs_length']
    sampler_config['baseline_tags'] = ['baseline-default']
    sampler_config['baseline_priors'  ]['baseline-default'] = [ [ 1.0 for j in range(length) ] for i in range(5) ]
    sampler_config['baseline_weights' ]['baseline-default'] = 1.0

def parse_sampler_config(config_file, sampler_config):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)
    if not config_parser.has_section('TFBS-Sampler'):
        raise IOError("Invalid configuration file.")
    sampler_config['seq_file'] = config_parser.get('TFBS-Sampler', 'sequences')
    if config_parser.has_option('TFBS-Sampler', 'alpha'):
        sampler_config['alpha'] = float(config_parser.get('TFBS-Sampler', 'alpha'))
    if config_parser.has_option('TFBS-Sampler', 'discount'):
        sampler_config['discount'] = float(config_parser.get('TFBS-Sampler', 'discount'))
    if config_parser.has_option('TFBS-Sampler', 'lambda'):
        sampler_config['lambda'] = float(config_parser.get('TFBS-Sampler', 'lambda'))
    if config_parser.has_option('TFBS-Sampler', 'process-prior'):
        sampler_config['process_prior'] = config_parser.get('TFBS-Sampler', 'process-prior').strip()
    if config_parser.has_option('TFBS-Sampler', 'background-model'):
        sampler_config['background_model'] = config_parser.get('TFBS-Sampler', 'background-model').strip()
    if config_parser.has_option('TFBS-Sampler', 'background-alpha'):
        sampler_config['background_alpha'] = read_matrix(config_parser, 'TFBS-Sampler', 'background-alpha', float)
    if config_parser.has_option('TFBS-Sampler', 'background-context'):
        sampler_config['background_context'] = int(config_parser.get('TFBS-Sampler', 'background-context'))
    if config_parser.has_option('TFBS-Sampler', 'background-weights'):
        sampler_config['background_weights'] = config_parser.get('TFBS-Sampler', 'background-weights')
    if config_parser.has_option('TFBS-Sampler', 'tfbs-length'):
        sampler_config['tfbs_length'] = int(config_parser.get('TFBS-Sampler', 'tfbs-length'))
    if config_parser.has_option('TFBS-Sampler', 'construct-graph'):
        sampler_config['construct_graph'] = str2bool(config_parser.get('TFBS-Sampler', 'construct_graph'))
    if config_parser.has_option('TFBS-Sampler', 'population-size'):
        sampler_config['population_size'] = int(config_parser.get('TFBS-Sampler', 'population-size'))
    if config_parser.has_option('TFBS-Sampler', 'save'):
        sampler_config['save'] = config_parser.get('TFBS-Sampler', 'save')
    if config_parser.has_option('TFBS-Sampler', 'samples'):
        sampler_config['samples'] = tuple(map(int, config_parser.get('TFBS-Sampler', 'samples').split(":")))
    if config_parser.has_option('TFBS-Sampler', 'baseline-priors'):
        sampler_config['baseline_priors']  = {}
        sampler_config['baseline_weights'] = {}
        sampler_config['baseline_tags']   = read_vector(config_parser, 'TFBS-Sampler', 'baseline-priors', str)
        for prior_name in sampler_config['baseline_tags']:
            sampler_config['baseline_priors'][prior_name] = read_matrix(config_parser, 'TFBS-Sampler', prior_name, float)
            if config_parser.has_option('TFBS-Sampler', '%s_weight' % prior_name):
                sampler_config['baseline_weights'][prior_name] = float(config_parser.get('TFBS-Sampler', '%s_weight' % prior_name))
            else:
                sampler_config['baseline_weights'][prior_name] = 1.0
        normalize_dict(sampler_config['baseline_weights'])
    else:
        generate_baseline(sampler_config)
    if config_parser.has_option('TFBS-Sampler', 'socket-file'):
        sampler_config['socket_file'] = config_parser.get('TFBS-Sampler', 'socket-file').strip()
