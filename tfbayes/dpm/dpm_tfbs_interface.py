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

# allow to send ctrl-c to the library
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

# load interface
# ------------------------------------------------------------------------------

_lib = None

if   os.path.exists(os.path.dirname(__file__)+'/.libs/libtfbayes-dpm.so'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libtfbayes-dpm.so')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/libtfbayes-dpm.dylib'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libtfbayes-dpm.dylib')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/cygdpm-0.dll'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/cygtfbates-dpm-0.dll')
else:
     for libname in ['libtfbayes-dpm.so.0', 'cygtfbayes-dpm-0.dll', 'libtfbayes-dpm.0.dylib']:
          if not _lib:
               try:
                    _lib = cdll.LoadLibrary(libname)
               except: pass

if not _lib:
     raise OSError('Couldn\'t find dpm library.')

# structures
# ------------------------------------------------------------------------------

class VECTOR(Structure):
     _fields_ = [("size", c_ulong),
                 ("vec",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_ulong),
                 ("columns", c_ulong),
                 ("mat",     POINTER(POINTER(c_double)))]

class PARTITION(Structure):
     _fields_ = []

class OPTIONS(Structure):
     _fields_ = [("tfbs_length",         c_ulong),
                 ("alpha",               c_double),
                 ("discount",            c_double),
                 ("lambda_",             c_double),
                 ("construct_graph",     c_bool),
                 ("process_prior",       c_char_p),
                 ("background_model",    c_char_p),
                 ("background_alpha",    POINTER(MATRIX)),
                 ("background_context",  c_ulong),
                 ("background_weights",  c_char_p),
                 ("baseline_weights",    POINTER(VECTOR)),
                 ("baseline_priors",     POINTER(POINTER(MATRIX))),
                 ("baseline_tags",       POINTER(c_char_p)),
                 ("baseline_n",          c_ulong),
                 ("partition",           POINTER(PARTITION)),
                 ("population_size",     c_ulong),
                 ("socket_file",         c_char_p)]

# function prototypes
# ------------------------------------------------------------------------------

_lib._alloc_vector.restype  = POINTER(VECTOR)
_lib._alloc_vector.argtypes = [c_ulong]

_lib._alloc_matrix.restype  = POINTER(MATRIX)
_lib._alloc_matrix.argtypes = [c_ulong, c_ulong]

_lib._free_vector.restype   = None
_lib._free_vector.argtypes  = [POINTER(VECTOR)]

_lib._free_matrix.restype   = None
_lib._free_matrix.argtypes  = [POINTER(MATRIX)]

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
_lib._dpm_tfbs_sample.argtypes  = [c_uint, c_uint]

_lib._dpm_tfbs_optimize.restype   = None
_lib._dpm_tfbs_optimize.argtypes  = []

_lib._dpm_tfbs_save.restype     = None
_lib._dpm_tfbs_save.argtypes    = [c_char_p]

_lib._dpm_tfbs_free.restype     = None
_lib._dpm_tfbs_free.argtypes    = []

_lib._dpm_tfbs_get_posterior.restype  = POINTER(MATRIX)
_lib._dpm_tfbs_get_posterior.argtypes = []

_lib._dpm_partition_new.restype  = POINTER(PARTITION)
_lib._dpm_partition_new.argtypes = []

_lib._dpm_partition_add_component.restype  = None
_lib._dpm_partition_add_component.argtypes = [POINTER(PARTITION), c_char_p]

_lib._dpm_partition_add_index.restype  = None
_lib._dpm_partition_add_index.argtypes = [POINTER(PARTITION), c_int, c_int]

_lib._dpm_partition_free.restype  = None
_lib._dpm_partition_free.argtypes = [POINTER(PARTITION)]

# convert datatypes
# ------------------------------------------------------------------------------

def copy_vector_to_c(v, c_v):
     for i in range(0, c_v.contents.size):
          c_v.contents.vec[i] = v[i]

def copy_matrix_to_c(m, c_m):
     for i in range(0, c_m.contents.rows):
          for j in range(0, c_m.contents.columns):
               c_m.contents.mat[i][j] = m[i][j]

def get_vector(c_v):
     v = []
     for i in range(0, c_v.contents.size):
          v.append(c_v.contents.vec[i])
     return v

def get_matrix(c_m):
     m = []
     for i in range(0, c_m.contents.rows):
          m.append([])
          for j in range(0, c_m.contents.columns):
               m[i].append(c_m.contents.mat[i][j])
     return m

def generate_c_partition(partition):
     c_partition = _lib._dpm_partition_new()
     for dpm_subset in partition:
          _lib._dpm_partition_add_component(c_partition, dpm_subset.identifier)
          for index in dpm_subset.subset:
               _lib._dpm_partition_add_index(c_partition, index[0], index[1])
     return c_partition

# functions that interface with the library
# ------------------------------------------------------------------------------

def dpm_init(options, input_file, partition=None):
     # do some sanity checks
     if not len(options['baseline_priors']) == len(options['baseline_weights']):
          raise IOError('Length mismatch between baseline priors and weights.')
     if not len(options['baseline_priors']) == len(options['baseline_tags']):
          raise IOError('Length mismatch between baseline priors and names.')
     if len(options['baseline_priors']) == 0:
          raise IOError('Length baseline priors specified.')

     # initialize simple c_options
     c_options                              = _lib._dpm_tfbs_options()
     c_options.contents.alpha               = options['alpha']
     c_options.contents.discount            = options['discount']
     c_options.contents.lambda_             = options['lambda']
     c_options.contents.construct_graph     = options['construct_graph']
     c_options.contents.process_prior       = options['process_prior']
     c_options.contents.background_model    = options['background_model']
     c_options.contents.background_context  = options['background_context']
     c_options.contents.background_weights  = options['background_weights']
     c_options.contents.tfbs_length         = options['tfbs_length']
     c_options.contents.population_size     = options['population_size']
     c_options.contents.socket_file         = options['socket_file']
     c_options.contents.partition           = generate_c_partition(partition) if partition else None

     # copy background alpha pseudo counts
     c_options.contents.background_alpha    = _lib._alloc_matrix(len(options['background_alpha'][0]), len(options['background_alpha']))
     copy_matrix_to_c(map(list, zip(*options['background_alpha'])), c_options.contents.background_alpha)

     baseline_size                          = len(options['baseline_priors'])
     c_options.contents.baseline_n          = c_ulong(baseline_size)
     # copy baseline priors and weights
     c_baseline_priors                      = (baseline_size*POINTER(MATRIX))()
     c_options.contents.baseline_weights    = _lib._alloc_vector(baseline_size)
     for idx, (name, prior) in enumerate(options['baseline_priors'].iteritems()):
          # first the prior
          c_baseline_priors[idx] = _lib._alloc_matrix(len(prior[0]), len(prior))
          copy_matrix_to_c(map(list, zip(*prior)), c_baseline_priors[idx])
          # and now its weight
          c_options.contents.baseline_weights.contents.vec[idx] = options['baseline_weights'][name]
     c_options.contents.baseline_priors     = cast(c_baseline_priors, POINTER(POINTER(MATRIX)))
     # copy baseline tags
     c_baseline_tags                        = (baseline_size*c_char_p)()
     for idx, name in enumerate(options['baseline_tags']):
          c_baseline_tags[idx] = c_char_p(name)
     c_options.contents.baseline_tags       = cast(c_baseline_tags, POINTER(c_char_p))

     # call the library
     _lib._dpm_tfbs_init(input_file)

     # free everything
     _lib._free_vector(c_options.contents.baseline_weights)
     for i in range(baseline_size):
          _lib._free_matrix(c_baseline_priors[i])
     if c_options.contents.partition:
          _lib._dpm_partition_free(c_options.contents.partition)

def dpm_print():
     _lib._dpm_tfbs_print()

def dpm_num_clusters():
     return _lib._dpm_tfbs_num_clusters()

def dpm_cluster_assignments():
     result  = _lib._dpm_tfbs_cluster_assignments()
     cluster = get_matrix(result)
     _lib._free_matrix(result)
     return cluster

def dpm_hist_likelihood():
     result     = _lib._dpm_tfbs_hist_likelihood()
     likelihood = get_vector(result)
     _lib._free_vector(result)
     return likelihood

def dpm_hist_switches():
     result   = _lib._dpm_tfbs_hist_switches()
     switches = get_vector(result)
     _lib._free_vector(result)
     return switches

def dpm_sample(n, burnin):
     c_n      = c_uint(n)
     c_burnin = c_uint(burnin)
     _lib._dpm_tfbs_sample(c_n, c_burnin)

def dpm_optimize():
     _lib._dpm_tfbs_optimize()

def dpm_save(filename):
     c_filename = c_char_p(filename)
     _lib._dpm_tfbs_save(c_filename)

def dpm_free():
     _lib._dpm_tfbs_free()

def dpm_get_posterior():
     return get_matrix(_lib._dpm_tfbs_get_posterior())
