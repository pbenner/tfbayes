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

import sys
import getopt
import os
import ConfigParser
import numpy as np
import math

# convert to fasta format
# ------------------------------------------------------------------------------

def fimport_description(description):
    result = '>'
    length = len(description)
    for i in range(0, length):
        if i+1 == length:
            result = result + description[i]
        else:
            result = result + description[i] + '|'
    return result+'\n'

def fimport_sequence(sequence):
    block_length = 80
    return ''.join([ '%.*s\n' % (block_length, sequence[i*block_length:]) for i in range(0, (len(sequence)-1)/block_length + 1) ])

def fimport(description, sequence):
    fdescription = fimport_description(description)
    fsequence    = fimport_sequence(sequence)
    return fdescription + fsequence

# parse fasta format
# ------------------------------------------------------------------------------

class parser():
    def __init__(self, file_name):
        self.fp = open(file_name)
        while 1:
            line = self.fp.readline()
            if not line or line[0] == '>':
                self.prev_line = line
                break
            if not line.strip() == '':
                raise IOError('Invalid fasta format.')
    def __del__(self):
        self.fp.close()
    def read_sequence(self):
        if not self.prev_line:
            return [], ''
        if not self.prev_line[0] == '>':
            raise IOError('Invalid fasta format.')
        description = self.prev_line[1:].split('|')
        description = map(lambda s: s.strip(), description)
        sequence = ''
        while 1:
            line = self.fp.readline()
            if not line or line[0] == '>':
                self.prev_line = line
                break
            if not line.strip() == '':
                sequence = sequence + line.strip()
        return description, sequence
