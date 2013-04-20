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

from ..interface import *
from ..uipac     import *

from ..          import phylotree

# load interface
# ------------------------------------------------------------------------------

_lib = load_library('tfbayes-alignment', 0)

# structures
# ------------------------------------------------------------------------------

class ALIGNMENT(Structure):
     def __repr__(self):
          return "ALIGNMENT()"

# function prototypes
# ------------------------------------------------------------------------------

_lib.alignment_new.restype  = POINTER(ALIGNMENT)
_lib.alignment_new.argtypes = [c_ulong, POINTER(phylotree.PT_ROOT)]

_lib.alignment_set.restype  = None
_lib.alignment_set.argtypes = [POINTER(ALIGNMENT), c_char_p, POINTER(VECTOR)]

_lib.alignment_marginal_likelihood.restype  = POINTER(VECTOR)
_lib.alignment_marginal_likelihood.argtypes = [POINTER(phylotree.PT_ROOT), POINTER(ALIGNMENT), POINTER(VECTOR)]

_lib.alignment_scan.restype  = POINTER(VECTOR)
_lib.alignment_scan.argtypes = [POINTER(phylotree.PT_ROOT), POINTER(ALIGNMENT), POINTER(MATRIX)]

_lib.alignment_free.restype  = None
_lib.alignment_free.argtypes = [POINTER(ALIGNMENT)]

# alignment data type
# ------------------------------------------------------------------------------

class alignment_t():
     def __repr__(self):
          return "alignment_t()"
     # this alignment is initialized with an alignment from
     # Bio.AlignIO
     def __init__(self, alignment, tree):
          c_record = _lib._alloc_vector(alignment.get_alignment_length())
          self.c_alignment = _lib.alignment_new(alignment.get_alignment_length(), tree)
          self.c_tree = tree
          for record in alignment:
               c_taxon = c_char_p(record.id)
               copy_vector_to_c(DNA.code(record), c_record)
               _lib.alignment_set(self.c_alignment, c_taxon, c_record)
          _lib._free_vector(c_record)
     def __del__(self):
          _lib.alignment_free(self.c_alignment)
     def marginal_likelihood(self, prior):
          c_prior  = _lib._alloc_vector(len(prior))
          copy_vector_to_c(prior, c_prior)
          c_result = _lib.alignment_marginal_likelihood(self.c_tree, self.c_alignment, c_prior)
          result   = get_vector(c_result)
          _lib._free_vector(c_result)
          _lib._free_vector(c_prior)
          return result
     def scan(self, counts):
          c_counts = _lib._alloc_matrix(len(counts), len(counts[0]))
          copy_matrix_to_c(counts, c_counts)
          c_result = _lib.alignment_scan(self.c_tree, self.c_alignment, c_counts)
          result   = get_vector(c_result)
          _lib._free_vector(c_result)
          _lib._free_matrix(c_counts)
          return result

