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
     _fields_ = [("size", c_ulong),
                 ("vec",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_ulong),
                 ("columns", c_ulong),
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

class CXX_VECTOR(Structure):
     length = 0
     def __new__(cdx, arg1):
          _lib._cxx_vector_alloc.restype  = POINTER(CXX_VECTOR)
          _lib._cxx_vector_alloc.argtypes = [c_ulong]
          if   type(arg1) == int or type(arg1) == long:
               return _lib._cxx_vector_alloc(arg1).contents
          elif type(arg1) == list:
               return _lib._cxx_vector_alloc(len(arg1)).contents
          else:
               raise ValueError("Invalid argument")
     def __init__(self, arg1):
          if type(arg1) == list:
               for (i,d) in enumerate(arg1):
                    self[i] = d
     def __del__(self):
          _lib._cxx_vector_free.restype  = None
          _lib._cxx_vector_free.argtypes = [POINTER(CXX_VECTOR)]
          _lib._cxx_vector_free(self)
     def __getitem__(self, i):
          _lib._cxx_vector_read.restype  = c_double
          _lib._cxx_vector_read.argtypes = [POINTER(CXX_VECTOR), c_ulong]
          return _lib._cxx_vector_read(self, i)
     def __setitem__(self, i, d):
          _lib._cxx_vector_write.restype  = None
          _lib._cxx_vector_write.argtypes = [POINTER(CXX_VECTOR), c_ulong, c_double]
          return _lib._cxx_vector_write(self, i, d)
     def export(self):
          return [ self[i] for i in range(self.length) ]

class CXX_MATRIX(Structure):
     def __new__(cdx, arg1, arg2 = None):
          _lib._cxx_matrix_alloc.restype  = POINTER(CXX_MATRIX)
          _lib._cxx_matrix_alloc.argtypes = [c_ulong, c_ulong]
          if arg2 is None:
               matrix = arg1
               if len(matrix) == 0:
                    return _lib._cxx_matrix_alloc(0, 0).contents
               rows    = len(matrix)
               columns = len(matrix[0])
               return _lib._cxx_matrix_alloc(rows, columns).contents
          else:
               return _lib._cxx_matrix_alloc(arg1, arg2).contents
     def __init__(self, arg1, arg2 = None):
          if arg2 is None:
               for (i, column) in enumerate(arg1):
                    for (j, d) in enumerate(column):
                         self[i,j] = d
     def __del__(self):
          _lib._cxx_matrix_free.restype  = None
          _lib._cxx_matrix_free.argtypes = [POINTER(CXX_MATRIX)]
          _lib._cxx_matrix_free(self)
     def __getitem__(self, (i, j)):
          _lib._cxx_matrix_read.restype  = c_double
          _lib._cxx_matrix_read.argtypes = [POINTER(CXX_MATRIX), c_ulong, c_ulong]
          return _lib._cxx_matrix_read(self, i, j)
     def __setitem__(self, (i, j), d):
          _lib._cxx_matrix_write.restype  = None
          _lib._cxx_matrix_write.argtypes = [POINTER(CXX_MATRIX), c_ulong, c_ulong, c_double]
          return _lib._cxx_matrix_write(self, i, j, d)

# c++ strings
# ------------------------------------------------------------------------------

class CXX_STRING(Structure):
     def __new__(cdx, str):
          _lib._cxx_string_alloc.restype  = POINTER(CXX_STRING)
          _lib._cxx_string_alloc.argtypes = [c_char_p]
          result = _lib._cxx_string_alloc(str).contents
          result.__init__(str)
          return result
     def __init__(self, str):
          pass
     def __del__(self):
          _lib._cxx_string_free.restype  = None
          _lib._cxx_string_free.argtypes = [POINTER(CXX_STRING)]
          _lib._cxx_string_free(self)
