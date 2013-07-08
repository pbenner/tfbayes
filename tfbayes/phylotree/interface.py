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

_lib = load_library('tfbayes-phylotree', 0)

# structures
# ------------------------------------------------------------------------------

class PT_ROOT(Structure):
    def __repr__(self):
        return "PT_ROOT()"
    def __str__(self):
        # TODO: this creates a memory leak
        return _lib.pt_print(self)

class ALIGNMENT(Structure):
     def __repr__(self):
          return "ALIGNMENT()"

# function prototypes
# ------------------------------------------------------------------------------

_lib.pt_parse_file.restype  = POINTER(PT_ROOT)
_lib.pt_parse_file.argtypes = [c_char_p]

_lib.pt_clone.restype  = POINTER(PT_ROOT)
_lib.pt_clone.argtypes = [POINTER(PT_ROOT)]

_lib.pt_destroy.restype  = None
_lib.pt_destroy.argtypes = [POINTER(PT_ROOT)]

_lib.pt_print.restype  = c_char_p
_lib.pt_print.argtypes = [POINTER(PT_ROOT)]

_lib.pt_index.restype  = c_long
_lib.pt_index.argtypes = [POINTER(PT_ROOT), c_char_p]

_lib.pt_num_leaves.restype  = c_ulong
_lib.pt_num_leaves.argtypes = [POINTER(PT_ROOT)]

_lib.pt_leaf_name.restype  = c_char_p
_lib.pt_leaf_name.argtypes = [POINTER(PT_ROOT), c_ulong]

_lib.pt_expectation.restype  = POINTER(VECTOR)
_lib.pt_expectation.argtypes = [POINTER(PT_ROOT), POINTER(VECTOR), POINTER(VECTOR)]

_lib.pt_approximate.restype  = POINTER(VECTOR)
_lib.pt_approximate.argtypes = [POINTER(PT_ROOT), POINTER(VECTOR)]

_lib.pt_dkl_optimize.restype  = POINTER(VECTOR)
_lib.pt_dkl_optimize.argtypes = [POINTER(PT_ROOT), POINTER(VECTOR)]

#
# ------------------------------------------------------------------------------

def pt_parse_file(filename):
    c_filename = c_char_p(filename)
    pt_root = _lib.pt_parse_file(c_filename)
    return pt_root

def pt_clone(pt_node):
    return _lib.pt_clone(pt_node)

def pt_destroy(pt_node):
    _lib.pt_destroy(pt_node)

def pt_print(pt_root):
    print pt_root.contents

def pt_index(pt_root, name):
    c_name = c_char_p(name)
    result = _lib.pt_index(pt_root, c_name)
    if (result == -1):
        return None
    else:
        return result

def pt_num_leaves(pt_root):
    return _lib.pt_num_leaves(pt_root)

def pt_leaf_name(pt_root, leaf):
    return _lib.pt_leaf_name(pt_root, leaf)

def pt_create_map(pt_root):
    n_leaves  = pt_num_leaves(pt_root)
    leaf_map = {}
    for i in range(n_leaves):
        name = pt_leaf_name(pt_root, i)
        leaf_map[name] = i
    return leaf_map

def pt_expectation(pt_root, observations, prior):
     c_observations = _lib._alloc_vector(len(observations))
     c_prior        = _lib._alloc_vector(len(prior))
     copy_vector_to_c(observations, c_observations)
     copy_vector_to_c(prior, c_prior)
     c_result = _lib.pt_expectation(pt_root, c_observations, c_prior)
     result   = get_vector(c_result)
     _lib._free_vector(c_result)
     _lib._free_vector(c_prior)
     _lib._free_vector(c_observations)
     return result

def pt_approximate(pt_root, observations):
     c_observations = _lib._alloc_vector(len(observations))
     copy_vector_to_c(observations, c_observations)
     c_result = _lib.pt_approximate(pt_root, c_observations)
     result   = get_vector(c_result)
     _lib._free_vector(c_result)
     _lib._free_vector(c_observations)
     return result
