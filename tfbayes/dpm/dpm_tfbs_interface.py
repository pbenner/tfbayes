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

# utilities
# ------------------------------------------------------------------------------

def transpose(matrix):
     return map(list, zip(*matrix))

# partitions
# ------------------------------------------------------------------------------

class PARTITION(c_void_p):
     def __init__(self, partition = None):
          _lib._dpm_partition_new.restype  = c_void_p
          _lib._dpm_partition_new.argtypes = []
          self.value = _lib._dpm_partition_new()
          if partition == None:
               return
          for dpm_subset in partition:
               self.add_component(dpm_subset.identifier)
               for index in dpm_subset.subset:
                    self.add_index(index[0], index[1])
     def add_component(self, subset_tag):
          _lib._dpm_partition_add_component.restype  = None
          _lib._dpm_partition_add_component.argtypes = [PARTITION, c_char_p]
          _lib._dpm_partition_add_component(self, subset_tag)
     def add_index(self, sequence, position):
          _lib._dpm_partition_add_index.restype  = None
          _lib._dpm_partition_add_index.argtypes = [PARTITION, c_int, c_int]
          _lib._dpm_partition_add_index(self, sequence, position)
     def free(self):
          _lib._dpm_partition_free.restype  = None
          _lib._dpm_partition_free.argtypes = [PARTITION]
          _lib._dpm_partition_free(self)

class PARTITION_LIST(c_void_p):
     def __init__(self, partition_list = None):
          _lib._dpm_partition_list_new.restype  = c_void_p
          _lib._dpm_partition_list_new.argtypes = []
          self.value = _lib._dpm_partition_list_new()
          if partition_list == None:
               return
          for partition in partition_list:
               self.add_partition(PARTITION(partition))
     def add_partition(self, partition):
          _lib._dpm_partition_list_add_partition.restype  = None
          _lib._dpm_partition_list_add_partition.argtypes = [PARTITION_LIST, PARTITION]
          _lib._dpm_partition_list_add_partition(self, partition)
     def free(self):
          _lib._dpm_partition_list_free.restype  = None
          _lib._dpm_partition_list_free.argtypes = [PARTITION_LIST]
          _lib._dpm_partition_list_free(self)
     def mean(self):
          _lib._dpm_mean.restype  = c_size_t
          _lib._dpm_mean.argtypes = [PARTITION_LIST]
          return _lib._dpm_mean(self)
     def median(self):
          _lib._dpm_median.restype  = c_size_t
          _lib._dpm_median.argtypes = [PARTITION_LIST]
          return _lib._dpm_median(self)

# options
# ------------------------------------------------------------------------------

class OPTIONS(Structure):
     _fields_ = [("tfbs_length",         c_size_t),
                 ("alpha",               c_double),
                 ("discount",            c_double),
                 ("lambda_",             c_double),
                 ("initial_temperature", c_double),
                 ("process_prior",       CXX_STRING),
                 ("background_model",    CXX_STRING),
                 ("background_alpha",    CXX_MATRIX),
                 ("background_context",  c_size_t),
                 ("background_weights",  CXX_STRING),
                 ("baseline_weights",    CXX_VECTOR),
                 ("baseline_priors",     POINTER(CXX_MATRIX)),
                 ("baseline_tags",       POINTER(CXX_STRING)),
                 ("baseline_n",          c_size_t),
                 ("partition",           PARTITION),
                 ("population_size",     c_size_t),
                 ("socket_file",         CXX_STRING)]
     def __new__(cls, options, partition):
          return _lib._dpm_tfbs_options().contents
     def __init__(self, options, partition):
          self.alpha               = options['alpha']
          self.discount            = options['discount']
          self.lambda_             = options['lambda']
          self.initial_temperature = options['initial_temperature']
          self.process_prior       = CXX_STRING(options['process_prior'])
          self.background_model    = CXX_STRING(options['background_model'])
          self.background_context  = options['background_context']
          self.background_weights  = CXX_STRING(options['background_weights'])
          self.tfbs_length         = options['tfbs_length']
          self.population_size     = options['population_size']
          self.socket_file         = CXX_STRING(options['socket_file'])
          self.partition           = PARTITION(partition)
          self.background_alpha    = CXX_MATRIX(transpose(options['background_alpha']))
          # copy baseline
          self.baseline_n          = len(options['baseline_priors'])
          self.baseline_tags       = (self.baseline_n*CXX_STRING)()
          for idx, name in enumerate(options['baseline_tags']):
               self.baseline_tags[idx] = CXX_STRING(name)
          # copy baseline priors and weights
          self.baseline_priors     = (self.baseline_n*CXX_MATRIX)()
          self.baseline_weights    = CXX_VECTOR(self.baseline_n)
          for idx, (name, prior) in enumerate(options['baseline_priors'].iteritems()):
               # first the prior
               self.baseline_priors[idx] = CXX_MATRIX(transpose(prior))
               # and now its weight
               self.baseline_weights[idx] = options['baseline_weights'][name]

# function prototypes
# ------------------------------------------------------------------------------

_lib._dpm_tfbs_options.restype  = POINTER(OPTIONS)
_lib._dpm_tfbs_options.argtypes = []

_lib._dpm_tfbs_init.restype  = None
_lib._dpm_tfbs_init.argtypes = [c_char_p, c_char_p]

_lib._dpm_tfbs_num_clusters.restype  = c_size_t
_lib._dpm_tfbs_num_clusters.argtypes = []

_lib._dpm_tfbs_cluster_assignments.restype  = POINTER(MATRIX)
_lib._dpm_tfbs_cluster_assignments.argtypes = []

_lib._dpm_tfbs_hist_likelihood.restype  = CXX_MATRIX
_lib._dpm_tfbs_hist_likelihood.argtypes = []

_lib._dpm_tfbs_hist_switches.restype  = CXX_MATRIX
_lib._dpm_tfbs_hist_switches.argtypes = []

_lib._dpm_tfbs_print.restype  = None
_lib._dpm_tfbs_print.argtypes = []

_lib._dpm_tfbs_sample.restype  = None
_lib._dpm_tfbs_sample.argtypes = [c_size_t, c_size_t]

_lib._dpm_tfbs_optimize.restype  = None
_lib._dpm_tfbs_optimize.argtypes = []

_lib._dpm_tfbs_save.restype  = None
_lib._dpm_tfbs_save.argtypes = [c_char_p]

_lib._dpm_tfbs_free.restype  = None
_lib._dpm_tfbs_free.argtypes = []

# functions that interface with the library
# ------------------------------------------------------------------------------

def dpm_init(options, phylogenetic_input, alignment_input, partition=None):
     # do some sanity checks
     if not len(options['baseline_priors']) == len(options['baseline_weights']):
          raise IOError('Length mismatch between baseline priors and weights.')
     if not len(options['baseline_priors']) == len(options['baseline_tags']):
          raise IOError('Length mismatch between baseline priors and names.')
     if len(options['baseline_priors']) == 0:
          raise IOError('Length baseline priors specified.')

     # initialize simple c_options
     c_options = pointer(OPTIONS(options, partition))

     # call the library
     _lib._dpm_tfbs_init(phylogenetic_input, alignment_input)

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
     likelihood = result.export()
     result.free()
     return likelihood

def dpm_hist_switches():
     result   = _lib._dpm_tfbs_hist_switches()
     switches = result.export()
     result.free()
     return switches

def dpm_sample(n, burnin):
     _lib._dpm_tfbs_sample(n, burnin)

def dpm_optimize():
     _lib._dpm_tfbs_optimize()

def dpm_save(filename):
     _lib._dpm_tfbs_save(filename)

def dpm_free():
     _lib._dpm_tfbs_free()

def dpm_mean(partition_list):
     c_partition_list = PARTITION_LIST(partition_list)
     i = c_partition_list.mean()
     c_partition_list.free()
     return partition_list[i]

def dpm_median(partition_list):
     c_partition_list = PARTITION_LIST(partition_list)
     i = c_partition_list.median()
     c_partition_list.free()
     return partition_list[i]
