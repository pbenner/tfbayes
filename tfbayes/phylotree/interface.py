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

if   os.path.exists(os.path.dirname(__file__)+'/.libs/libphylotree.so'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libphylotree.so')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/libphylotree.dylib'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libphylotree.dylib')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/cygphylotree-0.dll'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/cygphylotree-0.dll')
else:
     for libname in ['libphylotree.so.0', 'cygphylotree-0.dll', 'libphylotree.0.dylib']:
          if not _lib:
               try:
                    _lib = cdll.LoadLibrary(libname)
               except: pass

if not _lib:
     raise OSError('Couldn\'t find phylotree library.')

# structures
# ------------------------------------------------------------------------------

class VECTOR(Structure):
     _fields_ = [("size", c_ulong),
                 ("vec",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_ulong),
                 ("columns", c_ulong),
                 ("mat",     POINTER(POINTER(c_double)))]

class PT_ROOT(Structure):
    def __repr__(self):
        return "PT_ROOT()"
    def __str__(self):
        # TODO: this creates a memory leak
        return _lib.pt_print(self)

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

_lib._free.restype  = None
_lib._free.argtypes = [POINTER(None)]

_lib.pt_parse_file.restype  = POINTER(PT_ROOT)
_lib.pt_parse_file.argtypes = [c_char_p]

_lib.pt_print.restype  = c_char_p
_lib.pt_print.argtypes = [POINTER(PT_ROOT)]

_lib.pt_index.restype  = c_long
_lib.pt_index.argtypes = [POINTER(PT_ROOT), c_char_p]

_lib.pt_num_leafs.restype  = c_ulong
_lib.pt_num_leafs.argtypes = [POINTER(PT_ROOT)]

_lib.pt_leaf_name.restype  = c_char_p
_lib.pt_leaf_name.argtypes = [POINTER(PT_ROOT), c_ulong]

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

#
# ------------------------------------------------------------------------------

def pt_parse_file(filename):
    c_filename = c_char_p(filename)
    pt_root = _lib.pt_parse_file(c_filename)
    return pt_root

def pt_print(pt_root):
    print pt_root.contents

def pt_index(pt_root, name):
    c_name = c_char_p(name)
    result = _lib.pt_index(pt_root, c_name)
    if (result == -1):
        return None
    else:
        return result

def pt_num_leafs(pt_root):
    return _lib.pt_num_leafs(pt_root)

def pt_leaf_name(pt_root, leaf):
    return _lib.pt_leaf_name(pt_root, leaf)

def pt_create_map(pt_root):
    n_leafs  = pt_num_leafs(pt_root)
    leaf_map = {}
    for i in range(n_leafs):
        name = pt_leaf_name(pt_root, i)
        leaf_map[name] = i
    return leaf_map

def pt_expectation(pt_root, observations, prior):
     c_observations = _lib._alloc_vector(len(observations))
     c_prior        = _lib._alloc_vector(len(prior))
     c_result       = _lib._alloc_vector(len(prior))
     copy_vector_to_c(observations, c_observations)
     copy_vector_to_c(prior, c_prior)
     _lib.pt_expectation(pt_root, c_observations, c_prior, c_result)
     result = get_vector(c_result)
     _lib._free_vector(c_result)
     return result
