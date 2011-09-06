#! /usr/bin/env python

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

class DNA:
    letters = 'GATC'
    @staticmethod
    def is_nucleotide(n):
        if (n == 'A' or n == 'a' or n == 'C' or n == 'c' or
            n == 'G' or n == 'g' or n == 'T' or n == 't'):
            return True
        else:
            return False
    @staticmethod
    def code(sequence):
        codebook = { 'G': 0, 'g': 0, 'A': 1, 'a': 1,
                     'T': 2, 't': 2, 'C': 3, 'c': 3 }
        if len(sequence) == 1:
            return codebook[sequence]
        else:
            return map(lambda x: codebook[x], sequence)
    @staticmethod
    def decode(numbers):
        codebook = { 0: 'G', 1: 'A', 2: 'T', 3: 'C' }
        if type(numbers) == int:
            return codebook[numbers]
        else:
            return ''.join(map(lambda x: codebook[x], numbers))
