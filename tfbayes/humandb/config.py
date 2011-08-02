#! /usr/bin/env python

# Copyright (C) 2010 Philipp Benner
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

import os.path
import numpy as np
import ConfigParser

def readVector(config, section, option, converter):
    vector_str = config.get(section, option)
    vector     = map(converter, vector_str.split(' '))
    return vector

def readMatrix(config, section, option, converter):
    matrix_str = config.get(section, option)
    matrix     = []
    for line in matrix_str.split('\n'):
        if line != '':
            matrix.append([converter(a) for a in line.split(' ')])
    return matrix

def readStrategy(config_parser, section, options):
    if config_parser.has_option(section, 'strategy'):
        options['strategy'] = config_parser.get(section, 'strategy')

def readAlpha(config_parser, events, bins, section, converter):
    alpha = []
    if config_parser.has_option(section, 'alpha'):
        alpha = readMatrix(config_parser, section, 'alpha', converter)
        if not len(alpha) == events or not len(alpha[0]) == bins:
            raise ValueError("Wrong number of alpha parameters.")
    else:
        alpha = np.repeat(1, events*bins).reshape(events,bins)
    return alpha

def readModelPrior(config_parser, n, section, converter):
    models = []
    mprior = list(np.repeat(1, n))
    if config_parser.has_option(section, 'mprior'):
        mprior = list(np.repeat(0, n))
        models = readVector(config_parser, section, 'mprior', converter)
        num_models = len(models)
        for model in models:
            mprior[model-1] = 1.0/num_models
    return mprior

def readScript(config_parser, section, dir):
    if config_parser.has_option(section, 'script'):
        file = config_parser.get(section, 'script')
        f = open(os.path.abspath(dir)+'/'+file, 'r')
        str = f.read()
        f.close()
        return str
    else:
        return None
