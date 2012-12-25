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

from ..interface import *

# load interface
# ------------------------------------------------------------------------------

_lib = interface.load_library('tfbayes-dpm', 0)

# function prototypes
# ------------------------------------------------------------------------------

_lib._dpm_gaussian_init.restype     = None
_lib._dpm_gaussian_init.argtypes    = [c_int, c_double, POINTER(MATRIX), POINTER(MATRIX), POINTER(VECTOR), POINTER(VECTOR)]

_lib._dpm_gaussian_num_clusters.restype  = c_uint
_lib._dpm_gaussian_num_clusters.argtypes = []

_lib._dpm_gaussian_means.restype  = POINTER(MATRIX)
_lib._dpm_gaussian_means.argtypes = []

_lib._dpm_gaussian_data.restype   = POINTER(MATRIX)
_lib._dpm_gaussian_data.argtypes  = []

_lib._dpm_gaussian_cluster_tags.restype  = POINTER(VECTOR)
_lib._dpm_gaussian_cluster_tags.argtypes = []

_lib._dpm_gaussian_cluster_elements.restype  = POINTER(MATRIX)
_lib._dpm_gaussian_cluster_elements.argtypes = [c_int]

_lib._dpm_gaussian_cluster_assignments.restype  = POINTER(VECTOR)
_lib._dpm_gaussian_cluster_assignments.argtypes = []

_lib._dpm_gaussian_hist_likelihood.restype  = POINTER(VECTOR)
_lib._dpm_gaussian_hist_likelihood.argtypes = []

_lib._dpm_gaussian_hist_switches.restype    = POINTER(VECTOR)
_lib._dpm_gaussian_hist_switches.argtypes   = []

_lib._dpm_gaussian_original_means.restype   = POINTER(MATRIX)
_lib._dpm_gaussian_original_means.argtypes  = []

_lib._dpm_gaussian_original_cluster_assignments.restype   = POINTER(VECTOR)
_lib._dpm_gaussian_original_cluster_assignments.argtypes  = []

_lib._dpm_gaussian_print.restype    = None
_lib._dpm_gaussian_print.argtypes   = []

_lib._dpm_gaussian_sample.restype   = None
_lib._dpm_gaussian_sample.argtypes  = [c_uint, c_uint]

_lib._dpm_gaussian_free.restype     = None
_lib._dpm_gaussian_free.argtypes    = []

_lib._dpm_gaussian_get_posterior.restype  = POINTER(MATRIX)
_lib._dpm_gaussian_get_posterior.argtypes = []

#
# ------------------------------------------------------------------------------

def dpm_init(n, alpha, cov, cov_0, mu_0, pi):
     c_cov   = _lib._alloc_matrix(len(cov), len(cov[0]))
     c_cov_0 = _lib._alloc_matrix(len(cov_0), len(cov_0[0]))
     c_mu_0  = _lib._alloc_vector(len(mu_0))
     c_pi    = _lib._alloc_vector(len(pi))
     copy_matrix_to_c(cov, c_cov)
     copy_matrix_to_c(cov_0, c_cov_0)
     copy_vector_to_c(mu_0, c_mu_0)
     copy_vector_to_c(pi, c_pi)
     _lib._dpm_gaussian_init(n, alpha, c_cov, c_cov_0, c_mu_0, c_pi)
     _lib._free_matrix(c_cov)
     _lib._free_matrix(c_cov_0)
     _lib._free_vector(c_mu_0)
     _lib._free_vector(c_pi)

def dpm_print():
     _lib._dpm_gaussian_print()

def dpm_num_clusters():
     return _lib._dpm_gaussian_num_clusters()

def dpm_cluster_assignments():
     result  = _lib._dpm_gaussian_cluster_assignments()
     cluster = map(int, get_vector(result))
     _lib._free_vector(result)
     return cluster

def dpm_cluster_tags():
     result = _lib._dpm_gaussian_cluster_tags()
     tags   = map(int, get_vector(result))
     _lib._free_vector(result)
     return tags

def dpm_cluster_elements(tag):
     result   = _lib._dpm_gaussian_cluster_elements(tag)
     elements = get_matrix(result)
     _lib._free_matrix(result)
     return elements

def dpm_hist_likelihood():
     result     = _lib._dpm_gaussian_hist_likelihood()
     likelihood = get_vector(result)
     _lib._free_matrix(result)
     return likelihood

def dpm_hist_switches():
     result   = _lib._dpm_gaussian_hist_switches()
     switches = get_vector(result)
     _lib._free_vector(result)
     return switches

def dpm_means():
     result = _lib._dpm_gaussian_means()
     means  = get_matrix(result)
     _lib._free_matrix(result)
     return means

def dpm_data():
     result = _lib._dpm_gaussian_data()
     data  = get_matrix(result)
     _lib._free_matrix(result)
     return data

def dpm_original_means():
     result = _lib._dpm_gaussian_original_means()
     means  = get_matrix(result)
     _lib._free_matrix(result)
     return means

def dpm_original_cluster_assignments():
     result = _lib._dpm_gaussian_original_cluster_assignments()
     tags   = map(int, get_vector(result))
     _lib._free_vector(result)
     return tags

def dpm_sample(n, burnin):
     _lib._dpm_gaussian_sample(n, burnin)

def dpm_free():
     _lib._dpm_gaussian_free()

def dpm_get_posterior():
     return get_matrix(_lib._dpm_gaussian_get_posterior())
