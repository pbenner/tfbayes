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

from ctypes import *

# load interface
# ------------------------------------------------------------------------------

_lib = None

if   os.path.exists(os.path.dirname(__file__)+'/.libs/libhumandb.so'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libhumandb.so')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/libhumandb.dylib'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libhumandb.dylib')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/cygbayesian-binning-0.dll'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/cygbayesian-binning-0.dll')
else:
     for libname in ['libhumandb.so.0', 'cyghumandb-0.dll', 'libhumandb.0.dylib']:
          if not _lib:
               try:
                    _lib = cdll.LoadLibrary(libname)
               except: pass

if not _lib:
     raise OSError('Couldn\'t find humandb library.')

# structures
# ------------------------------------------------------------------------------

class VECTOR(Structure):
     _fields_ = [("size", c_int),
                 ("vec",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_int),
                 ("columns", c_int),
                 ("mat",     POINTER(POINTER(c_double)))]

# function prototypes
# ------------------------------------------------------------------------------

_lib._allocVector.restype       = POINTER(VECTOR)
_lib._allocVector.argtypes      = [c_int]

_lib._allocMatrix.restype       = POINTER(MATRIX)
_lib._allocMatrix.argtypes      = [c_int, c_int]

_lib._freeVector.restype        = None
_lib._freeVector.argtypes       = [POINTER(VECTOR)]

_lib._freeMatrix.restype        = None
_lib._freeMatrix.argtypes       = [POINTER(MATRIX)]

_lib._free.restype              = None
_lib._free.argtypes             = [POINTER(None)]

_lib._hdb_init.restype          = None
_lib._hdb_init.argtypes         = [c_char_p]

_lib._hdb_free.restype          = None
_lib._hdb_free.argtypes         = []

_lib._hdb_create.restype        = POINTER(None)
_lib._hdb_create.argtypes       = [c_char_p]

_lib._hdb_open.restype          = POINTER(None)
_lib._hdb_open.argtypes         = [c_char_p]

_lib._hdb_open_ro.restype       = POINTER(None)
_lib._hdb_open_ro.argtypes      = [c_char_p]

_lib._hdb_close.restype         = None
_lib._hdb_close.argtypes        = [POINTER(None)]

_lib._hdb_load_maf.restype      = None
_lib._hdb_load_maf.argtypes     = [POINTER(None), c_char_p]

_lib._hdb_get_sequence.restype  = None
_lib._hdb_get_sequence.argtypes = [POINTER(None), c_int, c_int, c_char_p]

_lib._hdb_search.restype        = None
_lib._hdb_search.argtypes       = [POINTER(POINTER(None)), c_int, POINTER(c_char_p), c_char_p]

_lib._hdb_search_pwm.restype    = None
_lib._hdb_search_pwm.argtypes   = [POINTER(POINTER(None)), c_int, POINTER(c_char_p), POINTER(MATRIX), c_double]

_lib._hdb_count_codons.restype  = None
_lib._hdb_count_codons.argtypes = [POINTER(None), POINTER(c_long), c_long, POINTER(c_long)]

_lib._hdb_count_codons_upstream.restype  = None
_lib._hdb_count_codons_upstream.argtypes = [POINTER(None), c_long, POINTER(c_long)]

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

def hdb_init(program_name):
     c_program_name = c_char_p(program_name)
     _lib._hdb_init(c_program_name)

def hdb_free():
     _lib._hdb_free()

def hdb_create(db_file_name, db_name):
     c_db_file_name = c_char_p(db_file_name)
     c_db_name      = c_char_p(db_name)
     ret = _lib._hdb_create(c_db_file_name, c_db_name)
     if ret == None:
          raise IOError("Opening database failed.")
     return ret

def hdb_open(db_file_name, db_name):
     c_db_file_name = c_char_p(db_file_name)
     c_db_name      = c_char_p(db_name)
     ret = _lib._hdb_open(c_db_file_name, c_db_name)
     if ret == None:
          raise IOError("Opening database failed.")
     return ret

def hdb_open_ro(db_file_name, db_name):
     c_db_file_name = c_char_p(db_file_name)
     c_db_name      = c_char_p(db_name)
     ret = _lib._hdb_open_ro(c_db_file_name, c_db_name)
     if ret == None:
          raise IOError("Opening database failed.")
     return ret

def hdb_close(dbp):
     _lib._hdb_close(dbp)

def hdb_load_maf(dbp, maf):
     c_maf = c_char_p(maf)
     _lib._hdb_load_maf(dbp, c_maf)

def hdb_get_sequence(dbp, pos, num):
     c_buf = create_string_buffer(num+1)
     c_pos = c_int(pos)
     c_num = c_int(num)
     _lib._hdb_get_sequence(dbp, c_pos, c_num, c_buf)
     return c_buf.value

def hdb_get_sequence_pure(dbp, pos, num):
     c_buf = create_string_buffer(num+1)
     c_pos = c_int(pos)
     c_num = c_int(num)
     _lib._hdb_get_sequence_pure(dbp, c_pos, c_num, c_buf)
     return c_buf.value

def hdb_search_pwm(dbp_list, db_names, pwm, threshold):
     dbp_list_n   = len(dbp_list)
     c_dbp_list_n = c_int(dbp_list_n)
     c_dbp_list   = (c_void_p*dbp_list_n)()
     c_db_names   = (c_char_p*dbp_list_n)()
     c_pwm        = _lib._allocMatrix(len(pwm), len(pwm[0]))
     copyMatrixToC(pwm,  c_pwm)
     c_threshold  = c_double(threshold)
     for i in range(0, dbp_list_n):
          c_dbp_list[i] = dbp_list[i]
          c_db_names[i] = c_char_p(db_names[i])
     _lib._hdb_search_pwm(c_dbp_list, c_dbp_list_n, c_db_names, c_pwm, c_threshold)

def hdb_search(dbp_list, db_names, sequence):
     dbp_list_n   = len(dbp_list)
     c_dbp_list_n = c_int(dbp_list_n)
     c_dbp_list   = (c_void_p*dbp_list_n)()
     c_db_names   = (c_char_p*dbp_list_n)()
     for i in range(0, dbp_list_n):
          c_dbp_list[i] = dbp_list[i]
          c_db_names[i] = c_char_p(db_names[i])
     c_sequence   = c_char_p(sequence)
     _lib._hdb_search(c_dbp_list, c_dbp_list_n, c_db_names, c_sequence)

def hdb_count_codons(dbp, positions):
     c_result    = (c_long*64)(0)
     c_positions = (c_long*len(positions))()
     c_len       = c_long(len(positions))
     for i in range(0, len(positions)):
          c_positions[i] = positions[i]
     _lib._hdb_count_codons(dbp, c_positions, c_len, c_result)
     result = []
     for i in c_result:
          result.append(i)
     return np.array(result)

def hdb_count_codons_upstream(dbp, length):
     c_result = (c_long*64)(0)
     c_length = c_long(length)
     _lib._hdb_count_codons_upstream(dbp, c_length, c_result)
     result = []
     for i in c_result:
          result.append(i)
     return np.array(result)