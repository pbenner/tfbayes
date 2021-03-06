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

import copy
import ConfigParser
import re

from tools import *

from tfbayes.interface import *
from tfbayes.dpm       import *

# default sampler config
# ------------------------------------------------------------------------------

def default_sampler_config():
    sampler_config = tfbs_options_t()
    sampler_config.alignment_file       = ""
    sampler_config.phylogenetic_file    = ""
    sampler_config.alpha                = 0.05
    sampler_config.discount             = 0.0
    sampler_config._lambda_             = 0.01
    sampler_config.block_samples        = False
    sampler_config.block_samples_period = 1
    sampler_config.metropolis_proposals = 4
    sampler_config.optimize             = False
    sampler_config.optimize_period      = 1
    sampler_config.initial_temperature  = 10.0
    sampler_config.process_prior        = 'pitman-yor process'
    sampler_config.background_model     = 'independence-dirichlet'
    sampler_config.background_alpha     = [[1, 1, 1, 1, 5]]
    sampler_config.background_beta      = [10.0, 10.0]
    sampler_config.background_gamma     = [5.0, 0.2]
    sampler_config.background_cache     = ''
    sampler_config.background_context   = 2
    sampler_config.background_weights   = []
    sampler_config.population_size      = 1
    sampler_config.baseline_lengths     = []
    sampler_config.baseline_names       = []
    sampler_config.baseline_priors      = []
    sampler_config.baseline_weights     = []
    sampler_config.socket_file          = ""
    sampler_config.samples              = (1000,100)
    sampler_config.threads              = 1
    sampler_config.save                 = ""
    sampler_config.verbose              = 0
    return sampler_config

# parse config
# ------------------------------------------------------------------------------

def tr(m):
    return map(list, zip(*m))

def normalize_weights(v):
    s = sum(v)
    for i in range(len(v)):
        v[i] = v[i]/s

def parse_sampler_config(config_file, sampler_config):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)
    if not config_parser.has_section('TFBS-Sampler'):
        raise IOError("Invalid sampler configuration file.")
    sampler_config.alignment_file    = config_parser.get('TFBS-Sampler', 'alignment-file')
    sampler_config.phylogenetic_file = config_parser.get('TFBS-Sampler', 'phylogenetic-file')
    if config_parser.has_option('TFBS-Sampler', 'alpha'):
        sampler_config.alpha = float(config_parser.get('TFBS-Sampler', 'alpha'))
    if config_parser.has_option('TFBS-Sampler', 'discount'):
        sampler_config.discount = float(config_parser.get('TFBS-Sampler', 'discount'))
    if config_parser.has_option('TFBS-Sampler', 'lambda'):
        sampler_config._lambda_ = float(config_parser.get('TFBS-Sampler', 'lambda'))
    if config_parser.has_option('TFBS-Sampler', 'block-samples'):
        sampler_config.block_samples = str2bool(config_parser.get('TFBS-Sampler', 'block-samples'))
    if config_parser.has_option('TFBS-Sampler', 'block-samples-period'):
        sampler_config.block_samples_period = int(config_parser.get('TFBS-Sampler', 'block-samples-period'))
        if not sampler_config.block_samples_period >= 1:
            raise IOError("Illegal block-samples-period specified.")
    if config_parser.has_option('TFBS-Sampler', 'metropolis-proposals'):
        if not int(config_parser.get('TFBS-Sampler', 'metropolis-proposals')) >= 0:
            raise IOError("Illegal number of metropolis proposals.")
        sampler_config.metropolis_proposals = int(config_parser.get('TFBS-Sampler', 'metropolis-proposals'))
    if config_parser.has_option('TFBS-Sampler', 'optimize'):
        sampler_config.optimize = str2bool(config_parser.get('TFBS-Sampler', 'optimize'))
    if config_parser.has_option('TFBS-Sampler', 'optimize-period'):
        sampler_config.optimize_period = int(config_parser.get('TFBS-Sampler', 'optimize-period'))
        if not sampler_config.optimize_period >= 1:
            raise IOError("Illegal optimize-period specified.")
    if config_parser.has_option('TFBS-Sampler', 'initial-temperature'):
        sampler_config.initial_temperature = float(config_parser.get('TFBS-Sampler', 'initial-temperature'))
    if config_parser.has_option('TFBS-Sampler', 'process-prior'):
        sampler_config.process_prior = config_parser.get('TFBS-Sampler', 'process-prior').strip()
    if config_parser.has_option('TFBS-Sampler', 'background-model'):
        sampler_config.background_model = config_parser.get('TFBS-Sampler', 'background-model').strip()
    if config_parser.has_option('TFBS-Sampler', 'background-alpha'):
        def convert_pseudocount(x):
            if x == '*'      : return -1
            if float(x) < 0.0: raise IOError("Illegal background pseudocounts")
            return float(x)
        sampler_config.background_alpha = tr(read_matrix(config_parser, 'TFBS-Sampler', 'background-alpha', convert_pseudocount))
    if config_parser.has_option('TFBS-Sampler', 'background-beta'):
        sampler_config.background_beta = read_vector(config_parser, 'TFBS-Sampler', 'background-beta', float)
        if (len(sampler_config.background_beta) != 2):
            raise IOError("Invalid background beta parameters")
    if config_parser.has_option('TFBS-Sampler', 'background-gamma'):
        sampler_config.background_gamma = read_vector(config_parser, 'TFBS-Sampler', 'background-gamma', float)
        if (len(sampler_config.background_gamma) != 2):
            raise IOError("Invalid background gamma parameters")
    if config_parser.has_option('TFBS-Sampler', 'background-cache'):
        sampler_config.background_cache = config_parser.get('TFBS-Sampler', 'background-cache').strip()
    if config_parser.has_option('TFBS-Sampler', 'background-context'):
        sampler_config.background_context = config_parser.get('TFBS-Sampler', 'background-context')
    if config_parser.has_option('TFBS-Sampler', 'background-weights'):
        sampler_config.background_weights = read_vector(config_parser, 'TFBS-Sampler', 'background-weights', float)
    if config_parser.has_option('TFBS-Sampler', 'median-partition'):
        sampler_config.median_partition = str2bool(config_parser.get('TFBS-Sampler', 'median-partition'))
    if config_parser.has_option('TFBS-Sampler', 'population-size'):
        sampler_config.population_size = int(config_parser.get('TFBS-Sampler', 'population-size'))
    if config_parser.has_option('TFBS-Sampler', 'threads'):
        sampler_config.threads = int(config_parser.get('TFBS-Sampler', 'threads'))
        if (sampler_config.threads <= 0):
            raise IOError("Invalid number of threads")
    if config_parser.has_option('TFBS-Sampler', 'save'):
        sampler_config.save = config_parser.get('TFBS-Sampler', 'save')
    if config_parser.has_option('TFBS-Sampler', 'samples'):
        sampler_config.samples = map(int, config_parser.get('TFBS-Sampler', 'samples').split(":"))
    if config_parser.has_option('TFBS-Sampler', 'baseline-priors'):
        sampler_config.baseline_priors  = []
        sampler_config.baseline_weights = []
        sampler_config.baseline_names   = read_vector(config_parser, 'TFBS-Sampler', 'baseline-priors', str)
        sampler_config.baseline_length  = []
        for prior_name in sampler_config.baseline_names:
            # read the pseudocounts matrix
            sampler_config.baseline_priors.append(tr(read_matrix(config_parser, 'TFBS-Sampler', prior_name, float)))
            # the baseline weights if they are specified
            if config_parser.has_option('TFBS-Sampler', '%s_weight' % prior_name):
                sampler_config.baseline_weights.append(config_parser.get('TFBS-Sampler', '%s_weight' % prior_name))
            else:
                sampler_config.baseline_weights.append(1.0)
            # and the baseline length
            if not config_parser.has_option('TFBS-Sampler', '%s_length' % prior_name):
                raise IOError("Baseline length required for `%s'" % prior_name)
            # the length can be given as a range
            m = re.match("\s*([0-9]+)\s*-\s*([0-9]+)\s*",
                         config_parser.get('TFBS-Sampler', '%s_length' % prior_name))
            # if length is a range
            if m:
                sampler_config.baseline_lengths.append(range(int(m.group(1)), int(m.group(2))+1))
            # otherwise it is a vector of lengths
            else:
                sampler_config.baseline_lengths.append(
                    read_vector(config_parser, 'TFBS-Sampler', '%s_length' % prior_name, int))
        normalize_weights(sampler_config.baseline_weights)
    else:
        generate_baseline(sampler_config)
    if config_parser.has_option('TFBS-Sampler', 'socket-file'):
        sampler_config.socket_file = config_parser.get('TFBS-Sampler', 'socket-file').strip()
    return sampler_config
