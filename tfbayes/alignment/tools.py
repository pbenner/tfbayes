#! /usr/bin/env python

# Copyright (C) 2012 Philipp Benner
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

import sys

def find_record(species, alignment):
    """Return the index of the record that is associated with the given species."""
    for idx, track in enumerate(alignment):
        if track.name == species:
            return idx
    sys.stderr.write('Warning, reference species not found.\n')
    return None

# to scan a sequence we need to select subsequences...
# ------------------------------------------------------------------------------

def select_subsequence(sequence, i, length, skip_gaps=False):
    if skip_gaps:
        # skip all subsequences that start with a gap
        if sequence[i] == '-':
            return None
        result = ""
        for j in range(i, len(sequence)):
            if not sequence[j] == '-':
                result += sequence[j]
            if len(result) == length:
                break
        # the resulting subsequence might be shorter
        # than what we need!
        if not len(result) is length:
            return None
        else:
            return result, j
    else:
        return sequence[i:i+length], i+length
