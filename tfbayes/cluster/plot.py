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

import numpy as np
import math

from ..corebio.seq         import Alphabet
from ..uipac.alphabet      import DNA
from ..weblogo             import *
from ..weblogo.colorscheme import nucleotide

import information

# cluster plotting
# ------------------------------------------------------------------------------

def plot_motif(motif, file_name, title, fout=None):
    print 'Generating %s ...' % file_name
    # convert probabilities to counts by simply scaling the pobabilities,
    # this is ok, since the logo is scaled by the discrete entropy and not
    # the differential entropy of the Dirichlet-compound posterior distribution
    counts   = [ [ int(round(motif[j][i]*1000)) for j in range(len(motif)) ] for i in range(len(motif[0])) ]
    alphabet = Alphabet(DNA.letters, zip(DNA.letters.lower(), DNA.letters))
    data     = LogoData.from_counts(alphabet, np.array(counts))
    options  = LogoOptions()
    options.color_scheme    = nucleotide
    options.logo_title      = title
    options.creator_text    = ''
    options.fineprint       = ''
    options.stacks_per_line = 60
    options.yaxis_scale     = 1.0
    options.scale_width     = False
    # use the natural logarithm as unit
    options.unit_name       = "nats"
    options.yaxis_label     = "p"
    # use discrete entropy to scale the logo
    # note: somehow the weblogo library requires that the entropy is inverted:
    # Log(|A|) - H(X)
    # that is, the stack height is simply defined by what is given by data.entropy,
    # since we use nats, the entropy should also be in the interval [0, 1]!
    data.entropy = map(lambda x: 1.0 - x/math.log(4), information.entropy(motif, len(motif), len(motif[0])))
    data.entropy_interval = None
    format = LogoFormat(data, options)
    if not fout:
        fout = open(file_name, 'w')
        pdf_formatter(data, format, fout)
        fout.close()
    else:
        pdf_formatter(data, format, fout)

def plot_counts(counts, alpha, file_name, title, fout=None):
    print 'Generating %s ...' % file_name
    # transpose the counts matrix and add the prior
    counts   = [ [ counts[j][i] + alpha[j][i] for j in range(len(counts)) ] for i in range(len(counts[0])) ]
    alphabet = Alphabet(DNA.letters, zip(DNA.letters.lower(), DNA.letters))
    # do not use the prior option from LogoData, it computes a strange mean relative entropy,
    # but by adding the pseudo counts to the counts we get what we want
    data     = LogoData.from_counts(alphabet, np.array(counts))
    options  = LogoOptions()
    options.color_scheme    = nucleotide
    options.logo_title      = title
    options.creator_text    = ''
    options.fineprint       = ''
    options.stacks_per_line = 60
    options.title_fontsize  = 7
    format = LogoFormat(data, options)
    if not fout:
        fout = open(file_name, 'w')
        pdf_formatter(data, format, fout)
        fout.close()
    else:
        pdf_formatter(data, format, fout)

def plot_cluster(cluster, basename, is_revcomp=False):
    # use a different file name and title if cluster is a reverse complement
    # of some other cluster
    if is_revcomp:
        title = 'Cluster %d (rc). A = %.1f. C = %d. RS = %.2f' % (cluster.identifier, cluster.average_counts(), cluster.components, cluster.r_sequence())
        file_name = '%s_cluster_%d_revcomp.pdf' % (basename, cluster.identifier)
    else:
        title = 'Cluster %d. A = %.1f. C = %d. RS = %.2f' % (cluster.identifier, cluster.average_counts(), cluster.components, cluster.r_sequence())
        file_name = '%s_cluster_%d.pdf' % (basename, cluster.identifier)
    # do the actual plotting here
    plot_counts(cluster.counts, cluster.alpha, file_name, title)
    # and return the file name of the pdf
    return file_name

def plot_cluster_list(cluster_list, basename, revcomp=False):
    files = []
    for cluster in cluster_list:
        files.append(plot_cluster(cluster, basename, False))
        if revcomp:
            files.append(plot_cluster(cluster.revcomp(), basename, True))
    return files
