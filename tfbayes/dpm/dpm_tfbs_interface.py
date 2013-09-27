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

from ..interface import *

# load library
# ------------------------------------------------------------------------------

_lib = load_library('tfbayes-dpm', 0)

# structures
# ------------------------------------------------------------------------------

class PARTITION(Structure):
     _fields_ = []

class PARTITION_LIST(Structure):
     _fields_ = []

class OPTIONS(Structure):
     _fields_ = [("tfbs_length",         c_ulong),
                 ("alpha",               c_double),
                 ("discount",            c_double),
                 ("lambda_",             c_double),
                 ("initial_temperature", c_double),
                 ("process_prior",       POINTER(CXX_STRING)),
                 ("background_model",    POINTER(CXX_STRING)),
                 ("background_alpha",    POINTER(CXX_MATRIX)),
                 ("background_context",  c_ulong),
                 ("background_weights",  POINTER(CXX_STRING)),
                 ("baseline_weights",    POINTER(VECTOR)),
                 ("baseline_priors",     POINTER(POINTER(MATRIX))),
                 ("baseline_tags",       POINTER(c_char_p)),
                 ("baseline_n",          c_ulong),
                 ("partition",           POINTER(PARTITION)),
                 ("population_size",     c_ulong),
                 ("socket_file",         POINTER(CXX_STRING))]
     def __new__(cls, options, partition):
          c_options =  _lib._dpm_tfbs_options().contents
          c_options.__init__(options, partition)
          return c_options
     def __init__(self, options, partition):
          self.alpha               = options['alpha']
          self.discount            = options['discount']
          self.lambda_             = options['lambda']
          self.initial_temperature = options['initial_temperature']
          self.process_prior       = pointer(CXX_STRING(options['process_prior']))
          self.background_model    = pointer(CXX_STRING(options['background_model']))
          self.background_context  = options['background_context']
          self.background_weights  = pointer(CXX_STRING(options['background_weights']))
          self.tfbs_length         = options['tfbs_length']
          self.population_size     = options['population_size']
          self.socket_file         = pointer(CXX_STRING(options['socket_file']))
          self.partition           = generate_c_partition(partition) if partition else None
          self.background_alpha    = pointer(CXX_MATRIX(map(list, zip(*options['background_alpha']))))

# function prototypes
# ------------------------------------------------------------------------------

_lib._dpm_tfbs_options.restype       = POINTER(OPTIONS)
_lib._dpm_tfbs_options.argtypes      = []

_lib._dpm_tfbs_init.restype          = None
_lib._dpm_tfbs_init.argtypes         = [c_char_p, c_char_p]

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

_lib._dpm_partition_list_new.restype  = POINTER(PARTITION_LIST)
_lib._dpm_partition_list_new.argtypes = []

_lib._dpm_partition_list_add_partition.restype  = None
_lib._dpm_partition_list_add_partition.argtypes = [POINTER(PARTITION_LIST), POINTER(PARTITION)]

_lib._dpm_partition_list_free.restype  = None
_lib._dpm_partition_list_free.argtypes = [POINTER(PARTITION_LIST)]

# global variables
# ------------------------------------------------------------------------------

c_options = None

# convert datatypes
# ------------------------------------------------------------------------------

def generate_c_partition(partition):
     c_partition = _lib._dpm_partition_new()
     for dpm_subset in partition:
          _lib._dpm_partition_add_component(c_partition, dpm_subset.identifier)
          for index in dpm_subset.subset:
               _lib._dpm_partition_add_index(c_partition, index[0], index[1])
     return c_partition

# functions that interface with the library
# ------------------------------------------------------------------------------

def dpm_init(options, phylogenetic_input, alignment_input, partition=None):
     global c_options
     # do some sanity checks
     if not len(options['baseline_priors']) == len(options['baseline_weights']):
          raise IOError('Length mismatch between baseline priors and weights.')
     if not len(options['baseline_priors']) == len(options['baseline_tags']):
          raise IOError('Length mismatch between baseline priors and names.')
     if len(options['baseline_priors']) == 0:
          raise IOError('Length baseline priors specified.')

     # initialize simple c_options
     c_options = pointer(OPTIONS(options, partition))

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
     _lib._dpm_tfbs_init(phylogenetic_input, alignment_input)

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
