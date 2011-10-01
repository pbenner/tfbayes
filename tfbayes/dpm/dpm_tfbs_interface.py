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

import os
import numpy as np
import math
import ConfigParser
import getopt

from ctypes import *

# load interface
# ------------------------------------------------------------------------------

_lib = None

if   os.path.exists(os.path.dirname(__file__)+'/.libs/libdpm.so'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libdpm.so')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/libdpm.dylib'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libdpm.dylib')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/cygdpm-0.dll'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/cygdpm-0.dll')
else:
     for libname in ['libdpm.so.0', 'cygdpm-0.dll', 'libdpm.0.dylib']:
          if not _lib:
               try:
                    _lib = cdll.LoadLibrary(libname)
               except: pass

if not _lib:
     raise OSError('Couldn\'t find dpm library.')

# structures
# ------------------------------------------------------------------------------

class VECTOR(Structure):
     _fields_ = [("size", c_int),
                 ("vec",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_int),
                 ("columns", c_int),
                 ("mat",     POINTER(POINTER(c_double)))]

class OPTIONS(Structure):
     _fields_ = [("tfbs_length",       c_int),
                 ("alpha",             c_double),
                 ("discount",          c_double),
                 ("lambda_",           c_double),
                 ("process_prior",     c_char_p),
                 ("population_size",   c_int),
                 ("baseline_weights",  POINTER(VECTOR)),
                 ("baseline_priors",   POINTER(POINTER(MATRIX))),
                 ("baseline_n",        c_int)]

# function prototypes
# ------------------------------------------------------------------------------

_lib._allocVector.restype  = POINTER(VECTOR)
_lib._allocVector.argtypes = [c_int]

_lib._allocMatrix.restype  = POINTER(MATRIX)
_lib._allocMatrix.argtypes = [c_int, c_int]

_lib._freeVector.restype   = None
_lib._freeVector.argtypes  = [POINTER(VECTOR)]

_lib._freeMatrix.restype   = None
_lib._freeMatrix.argtypes  = [POINTER(MATRIX)]

_lib._free.restype         = None
_lib._free.argtypes        = [POINTER(None)]

_lib._dpm_tfbs_options.restype       = POINTER(OPTIONS)
_lib._dpm_tfbs_options.argtypes      = []

_lib._dpm_tfbs_init.restype          = None
_lib._dpm_tfbs_init.argtypes         = [c_char_p]

_lib._dpm_tfbs_num_clusters.restype  = c_uint
_lib._dpm_tfbs_num_clusters.argtypes = []

_lib._dpm_tfbs_cluster_assignments.restype  = POINTER(MATRIX)
_lib._dpm_tfbs_cluster_assignments.argtypes = []

_lib._dpm_tfbs_hist_likelihood.restype  = POINTER(VECTOR)
_lib._dpm_tfbs_hist_likelihood.argtypes = []

_lib._dpm_tfbs_hist_switches.restype  = POINTER(VECTOR)
_lib._dpm_tfbs_hist_switches.argtypes = []

_lib._dpm_tfbs_print.restype    = None
_lib._dpm_tfbs_print.argtypes   = []

_lib._dpm_tfbs_sample.restype   = None
_lib._dpm_tfbs_sample.argtypes  = []

_lib._dpm_tfbs_save.restype     = None
_lib._dpm_tfbs_save.argtypes    = [c_char_p]

_lib._dpm_tfbs_free.restype     = None
_lib._dpm_tfbs_free.argtypes    = []

_lib._dpm_tfbs_get_posterior.restype  = POINTER(MATRIX)
_lib._dpm_tfbs_get_posterior.argtypes = []

# convert datatypes
# ------------------------------------------------------------------------------

def copyVectorToC(v, c_v):
     for i in range(0, c_v.contents.size):
          c_v.contents.vec[i] = v[i]

def copyMatrixToC(m, c_m):
     for i in range(0, c_m.contents.rows):
          for j in range(0, c_m.contents.columns):
               c_m.contents.mat[i][j] = m[i][j]

def getVector(c_v):
     v = []
     for i in range(0, c_v.contents.size):
          v.append(c_v.contents.vec[i])
     return v

def getMatrix(c_m):
     m = []
     for i in range(0, c_m.contents.rows):
          m.append([])
          for j in range(0, c_m.contents.columns):
               m[i].append(c_m.contents.mat[i][j])
     return m

#
# ------------------------------------------------------------------------------

def dpm_init(options, input_file):
     c_options = _lib._dpm_tfbs_options()
     c_options.contents.alpha = options['alpha']
     c_options.contents.discount = options['discount']
     c_options.contents.lambda_ = options['lambda']
     c_options.contents.process_prior = options['process_prior']
     c_options.contents.tfbs_length = options['tfbs_length']
     c_options.contents.population_size = options['population_size']
     c_input_file = c_char_p(input_file)
     if not len(options['baseline_priors']) == len(options['baseline_weights']):
          raise IOError('Length mismatch between baseline priors and weights')
     prior_length = len(options['baseline_priors'])
     if prior_length > 0:
          normalized_weights = map(lambda x: float(x)/sum(options['baseline_weights']), options['baseline_weights'])
          c_baseline_weights = _lib._allocVector(prior_length)
          copyVectorToC(normalized_weights, c_baseline_weights)
          c_baseline_priors = (prior_length*POINTER(MATRIX))()
          for i, prior in zip(range(prior_length), options['baseline_priors']):
               c_baseline_priors[i] = _lib._allocMatrix(len(prior[0]), len(prior))
               copyMatrixToC(map(list, zip(*prior)), c_baseline_priors[i])
          c_options.contents.baseline_weights = c_baseline_weights
          c_options.contents.baseline_priors = pointer(c_baseline_priors[0])
          c_options.contents.baseline_n = c_int(prior_length)
     _lib._dpm_tfbs_init(c_input_file)
     _lib._freeVector(c_baseline_weights)
     for i in range(prior_length):
          _lib._freeMatrix(c_baseline_priors[i])

def dpm_print():
     _lib._dpm_tfbs_print()

def dpm_num_clusters():
     return _lib._dpm_tfbs_num_clusters()

def dpm_cluster_assignments():
     result  = _lib._dpm_tfbs_cluster_assignments()
     cluster = getMatrix(result)
     _lib._freeMatrix(result)
     return cluster

def dpm_hist_likelihood():
     result     = _lib._dpm_tfbs_hist_likelihood()
     likelihood = getVector(result)
     _lib._freeVector(result)
     return likelihood

def dpm_hist_switches():
     result   = _lib._dpm_tfbs_hist_switches()
     switches = getVector(result)
     _lib._freeVector(result)
     return switches

def dpm_sample(n, burnin):
     c_n      = c_int(n)
     c_burnin = c_int(burnin)
     _lib._dpm_tfbs_sample(c_n, c_burnin)

def dpm_save(filename):
     c_filename = c_char_p(filename)
     _lib._dpm_tfbs_save(c_filename)

def dpm_free():
     _lib._dpm_tfbs_free()

def dpm_get_posterior():
     return getMatrix(_lib._dpm_tfbs_get_posterior())
