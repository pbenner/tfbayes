#! /usr/bin/env python

# Copyright (C) 2011, 2012 Philipp Benner
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

from ctypes    import *
from datatypes import *

def load_library(name, major_version):
    lib = None

    if   os.path.exists(os.path.dirname(__file__)+'/.libs/lib%s.so' % name):
        lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/lib%s.so' % name)
    elif os.path.exists(os.path.dirname(__file__)+'/.libs/lib%s.dylib' % name):
        lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/lib%s.dylib' % name)
    elif os.path.exists(os.path.dirname(__file__)+'/.libs/cyg%s-%d.dll' % (name, major_version)):
        lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/cyg%s-%d.dll' % (name, major_version))
    else:
        for filename in ['lib%s.so.%d'    % (name, major_version),
                         'cyg%s-%d.dll'   % (name, major_version),
                         'lib%s.%d.dylib' % (name, major_version)]:
            if not lib:
                try:
                    lib = cdll.LoadLibrary(filename)
                except: pass

    if not lib:
        raise OSError('Could not find %s library.' % name)

    # these functions exist in every tfbayes library
    lib._alloc_vector.restype  = POINTER(VECTOR)
    lib._alloc_vector.argtypes = [c_ulong]
    lib._alloc_matrix.restype  = POINTER(MATRIX)
    lib._alloc_matrix.argtypes = [c_ulong, c_ulong]
    lib._free_vector.restype   = None
    lib._free_vector.argtypes  = [POINTER(VECTOR)]
    lib._free_matrix.restype   = None
    lib._free_matrix.argtypes  = [POINTER(MATRIX)]
    lib._free.restype          = None
    lib._free.argtypes         = [POINTER(None)]

    return lib
