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

import Bio.Motif as Motif

from ..motif.information import r_sequence
from ..uipac.alphabet import DNA

# motif plotting
# ------------------------------------------------------------------------------

class MotifStream():
    def __init__(self, counts):
        self.counts = []
        for line in counts:
            self.counts.append(' '.join(map(str, line)))
        self.counts.append('\n')
        self.counts.append('\n')
        self.i      = 0
    def readline(self):
        if (self.i < len(self.counts)):
            line = self.counts[self.i]
            self.i += 1
            return line
        else:
            return None

def plot_motif(motif, file_name, title):
    print 'Generating %s...' % file_name
    counts = [ [ int(round(motif[i][j]*1000)) for j in range(0, len(motif[0])) ] for i in range(0, len(motif)) ]
    stream = MotifStream(counts)
    m = Motif.Motif()
    m._from_horiz_matrix(stream, letters=DNA.letters)
    m.weblogo(file_name, title=title)

def plot_motifs(motifs, components, basename):
    for n, motif, comp in zip(range(0, len(motifs)), motifs, components):
        r_seq = r_sequence(motif)
        file_name = '%s_cluster_%d.png' % (basename, n)
        title = 'cluster_%d:%d, R_seq = %f' % (n, comp, r_seq)
        plot_motif(motif, file_name, title)
