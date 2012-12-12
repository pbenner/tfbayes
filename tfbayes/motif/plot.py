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

import numpy as np
import math

from corebio.seq import Alphabet
from weblogolib import *
from weblogolib.colorscheme import nucleotide

from ..motif.tools import reverse_complement
from ..motif.information import r_sequence
from ..uipac.alphabet import DNA

# motif plotting
# ------------------------------------------------------------------------------

def plot_probability(motif, file_name, title, fout=None):
    print 'Generating %s...' % file_name
    counts = [ [ int(round(motif[j][i]*1000)) for j in range(len(motif)) ] for i in range(len(motif[0])) ]
    alphabet = Alphabet(DNA.letters, zip(DNA.letters.lower(), DNA.letters))
    data = LogoData.from_counts(alphabet, np.array(counts))
    # use entropy to scale the logo
    for i in range(len(data.entropy)):
        data.entropy[i] = sum([ motif[j][i] for j in range(len(motif)) ])
    options = LogoOptions()
    options.color_scheme = nucleotide
    options.logo_title = title
    options.creator_text = ''
    options.fineprint = ''
    options.stacks_per_line = 60
    options.yaxis_scale = 1.0
    options.scale_width = False
    options.unit_name = "nats"
    options.yaxis_label = "p"
    format = LogoFormat(data, options)
    if not fout:
        fout = open(file_name, 'w')
        pdf_formatter(data, format, fout)
        fout.close()
    else:
        pdf_formatter(data, format, fout)

def plot_motif(motif, file_name, title, fout=None):
    print 'Generating %s ...' % file_name
    counts = [ [ motif[j][i] for j in range(4) ] for i in range(len(motif[0])) ]
    alphabet = Alphabet(DNA.letters, zip(DNA.letters.lower(), DNA.letters))
    data = LogoData.from_counts(alphabet, np.array(counts))
    options = LogoOptions()
    options.color_scheme    = nucleotide
    options.logo_title      = title
    options.creator_text    = ''
    options.fineprint       = ''
    options.stacks_per_line = 60
    options.title_fontsize  = 8
    format = LogoFormat(data, options)
    if not fout:
        fout = open(file_name, 'w')
        pdf_formatter(data, format, fout)
        fout.close()
    else:
        pdf_formatter(data, format, fout)

def average_components(motif):
    counts  = [ sum([ motif[j][i] for j in range(4) ]) for i in range(len(motif[0])) ]
    average = sum(counts)/float(len(motif[0]))
    return average

def plot_motifs(motifs, components, basename, revcomp=False):
    files = []
    for n, motif, comp in zip(range(0, len(motifs)), motifs, components):
        file_name = '%s_cluster_%d.pdf' % (basename, n)
        title = 'Cluster %d. A = %.1f. C = %d' % (n, average_components(motif), comp)
        plot_motif(motif, file_name, title)
        files.append(file_name)
        if revcomp:
            motif_revcomp = reverse_complement(motif)
            file_name = '%s_cluster_%d_revcomp.pdf' % (basename, n)
            title = 'Cluster %d. A = %.1f. C = %d' % (n, average_components(motif), comp)
            plot_motif(motif_revcomp, file_name, title)
            files.append(file_name)
    return files
