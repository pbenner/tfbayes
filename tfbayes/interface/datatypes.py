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
from tools  import *

# load library
# ------------------------------------------------------------------------------

_lib = load_library('tfbayes-interface', 0)

# c++ vector/matrix class
# ------------------------------------------------------------------------------

class CXX_VECTOR_PTR(c_void_p):
     pass

_lib._cxx_vector_alloc.restype  = CXX_VECTOR_PTR
_lib._cxx_vector_alloc.argtypes = [c_ulong]

_lib._cxx_vector_free.restype  = None
_lib._cxx_vector_free.argtypes = [CXX_VECTOR_PTR]

_lib._cxx_vector_read.restype  = c_double
_lib._cxx_vector_read.argtypes = [CXX_VECTOR_PTR, c_ulong]

_lib._cxx_vector_write.restype  = None
_lib._cxx_vector_write.argtypes = [CXX_VECTOR_PTR, c_ulong, c_double]

class CXX_VECTOR(Structure):
     ptr    = None
     length = 0
     def __init__(self, arg1):
          if   type(arg1) == int:
               self.ptr    = _lib._cxx_vector_alloc(arg1)
               self.length = arg1
          elif type(arg1) == list:
               self.ptr    = _lib._cxx_vector_alloc(len(arg1))
               self.length = len(arg1)
               for (i,d) in enumerate(arg1):
                    self[i] = d
          else:
               raise ValueError("Invalid argument")
     def __del__(self):
          _lib._cxx_vector_free(self.ptr)
     def __getitem__(self, i):
          return _lib._cxx_vector_read(self.ptr, i)
     def __setitem__(self, i, d):
          return _lib._cxx_vector_write(self.ptr, i, d)
     def export(self):
          return [ self[i] for i in range(self.length) ]

class CXX_MATRIX_PTR(c_void_p):
     pass

_lib._cxx_matrix_alloc.restype  = CXX_MATRIX_PTR
_lib._cxx_matrix_alloc.argtypes = [c_ulong, c_ulong]

_lib._cxx_matrix_read.restype  = c_double
_lib._cxx_matrix_read.argtypes = [CXX_MATRIX_PTR, c_ulong, c_ulong]

_lib._cxx_matrix_write.restype  = None
_lib._cxx_matrix_write.argtypes = [CXX_MATRIX_PTR, c_ulong, c_ulong, c_double]

_lib._cxx_matrix_free.restype  = None
_lib._cxx_matrix_free.argtypes = [CXX_MATRIX_PTR]

class CXX_MATRIX(Structure):
     ptr     = None
     rows    = 0
     columns = 0
     def __init__(self, arg1, arg2 = None):
          if arg2 is None:
               matrix = arg1
               if len(matrix) == 0:
                    self.ptr = _lib._cxx_matrix_alloc(0, 0)
               else:
                    self.rows     = len(matrix)
                    self.columns  = len(matrix[0])
                    self.ptr = _lib._cxx_matrix_alloc(self.rows, self.columns)
                    # initialize matrix
                    for (i, column) in enumerate(matrix):
                         for (j, d) in enumerate(column):
                              self[i,j] = d
          else:
               self.rows    = arg1
               self.columns = arg2
               self.ptr     = _lib._cxx_matrix_alloc(self.rows, self.columns)
     def __del__(self):
          _lib._cxx_matrix_free(self.ptr)
     def __getitem__(self, (i, j)):
          return _lib._cxx_matrix_read(self.ptr, i, j)
     def __setitem__(self, (i, j), d):
          return _lib._cxx_matrix_write(self.ptr, i, j, d)
     def export(self):
          return [ [ self[i,j] for j in range(self.columns) ] for i in range(self.rows) ]

# c vector/matrix class
# ------------------------------------------------------------------------------

class VECTOR(Structure):
     _fields_ = [("size", c_ulong),
                 ("vec",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_ulong),
                 ("columns", c_ulong),
                 ("mat",     POINTER(POINTER(c_double)))]
