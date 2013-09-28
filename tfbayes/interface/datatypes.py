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

from ctypes import *

# c vector/matrix class
# ------------------------------------------------------------------------------

class VECTOR(Structure):
     _fields_ = [("size", c_size_t),
                 ("vec",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_size_t),
                 ("columns", c_size_t),
                 ("mat",     POINTER(POINTER(c_double)))]

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

# load library
# ------------------------------------------------------------------------------

from tfbayes.interface.loadlibrary import *

_lib = load_library('tfbayes-interface', 0)

# c++ vector/matrix class
# ------------------------------------------------------------------------------

class CXX_VECTOR(c_void_p):
     def __init__(self, arg1 = None):
          _lib._cxx_vector_alloc.restype  = c_void_p
          _lib._cxx_vector_alloc.argtypes = [c_size_t]
          if   type(arg1) == int or type(arg1) == long:
               self.value = _lib._cxx_vector_alloc(arg1)
          elif type(arg1) == list:
               self.value = _lib._cxx_vector_alloc(len(arg1))
               # init vector
               for (i,d) in enumerate(arg1):
                    self[i] = d
          else:
               raise ValueError("Invalid argument")
     def __getitem__(self, i):
          _lib._cxx_vector_read.restype  = c_double
          _lib._cxx_vector_read.argtypes = [c_void_p, c_size_t]
          return _lib._cxx_vector_read(self.value, i)
     def __setitem__(self, i, d):
          _lib._cxx_vector_write.restype  = None
          _lib._cxx_vector_write.argtypes = [c_void_p, c_size_t, c_double]
          return _lib._cxx_vector_write(self.value, i, d)
     def free(self):
          _lib._cxx_vector_free.restype  = None
          _lib._cxx_vector_free.argtypes = [c_void_p]
          _lib._cxx_vector_free(self.value)
     def export(self):
          return [ self[i] for i in range(self.length) ]

class CXX_MATRIX(c_void_p):
     def __init__(self, arg1, arg2 = None):
          _lib._cxx_matrix_alloc.restype  = c_void_p
          _lib._cxx_matrix_alloc.argtypes = [c_size_t, c_size_t]
          if arg2 is None:
               if len(arg1) == 0:
                    self.value = _lib._cxx_matrix_alloc(0, 0)
                    return
               rows    = len(arg1)
               columns = len(arg1[0])
               self.value = _lib._cxx_matrix_alloc(rows, columns)
               # initialize matrix
               for (i, column) in enumerate(arg1):
                    for (j, d) in enumerate(column):
                         self[i,j] = d
          else:
               self.value = _lib._cxx_matrix_alloc(arg1, arg2)
     def __getitem__(self, (i, j)):
          _lib._cxx_matrix_read.restype  = c_double
          _lib._cxx_matrix_read.argtypes = [self.value, c_size_t, c_size_t]
          return _lib._cxx_matrix_read(self.value, i, j)
     def __setitem__(self, (i, j), d):
          _lib._cxx_matrix_write.restype  = None
          _lib._cxx_matrix_write.argtypes = [c_void_p, c_size_t, c_size_t, c_double]
          return _lib._cxx_matrix_write(self.value, i, j, d)
     def free(self):
          _lib._cxx_matrix_free.restype  = None
          _lib._cxx_matrix_free.argtypes = [c_void_p]
          _lib._cxx_matrix_free(self.value)

# c++ strings
# ------------------------------------------------------------------------------

class CXX_STRING(c_void_p):
     def __init__(self, str):
          _lib._cxx_string_alloc.restype  = c_void_p
          _lib._cxx_string_alloc.argtypes = [c_char_p]
          self.value = _lib._cxx_string_alloc(str)
     def __str__(self):
          _lib._cxx_string_getstr.restype  = c_char_p
          _lib._cxx_string_getstr.argtypes = [c_void_p]
          result = _lib._cxx_string_getstr(self.value)
          return str(result)
     def free(self):
          _lib._cxx_string_free.restype  = None
          _lib._cxx_string_free.argtypes = [c_void_p]
          _lib._cxx_string_free(self.value)
