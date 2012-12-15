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

import re
import fasta

# parse sequences file (approximate multiple alignments)
# ------------------------------------------------------------------------------

def filter_line(line):
    line = re.sub(r'\s+', ' ', line).strip()
    line = line.split(';')
    line.remove('')
    line = map(lambda s: s.strip(), line)
    return line

def parse_line(lines):
    return map(lambda s:
                  map(lambda e: float(e), s.split(' ')),
              lines)

def parse_sequences(seq_file):
    parser = fasta.parser(seq_file)
    descriptions = []
    sequences    = []
    while 1:
        description, line = parser.read_sequence()
        if not line:
            break
        lines = filter_line(line)
        descriptions.append(description[0])
        sequences.append(parse_line(lines))
    return descriptions, sequences
