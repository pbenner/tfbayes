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

import math
import tools

# R_sequence
# ------------------------------------------------------------------------------

def r_sequence(cluster_counts):
    z = 1.0/float(len(cluster_counts[0]))
    cluster_freq = tools.compute_frequencies(cluster_counts)
    return z*sum([ 2 + sum([ cluster_freq[i][j]*math.log(cluster_freq[i][j], 2)
                             for i in range(0, len(cluster_freq)) ])
                   for j in range(0, len(cluster_freq[0])) ])

# information content
# ------------------------------------------------------------------------------

def information_content_column(j, bg_freq, motif):
    return sum([ motif[i][j]*math.log(motif[i][j]/bg_freq[i][0], 2) for i in range(0, len(motif)) ])

def information_content(bg_freq, motif):
    norm = 1.0/float(len(motif[0]))
    return norm*sum([ sum([ motif[i][j]*math.log(motif[i][j]/bg_freq[i][0], 2) for j in range(0, len(motif[0])) ]) for i in range(0, len(motif)) ])

# Kullback-Leibler divergence
# ------------------------------------------------------------------------------

def kl_divergence_column(j, k, bg_freq, c1_freq, c2_freq, length):
    j1 = j
    j2 = j-k
    if   j1 >= 0 and j1 < length and j2 >= 0 and j2 < length:
        # both matrices overlap
        return sum( [ c1_freq[i][j1]*math.log(c1_freq[i][j1]/c2_freq[i][j2], 2) for i in range(0, len(c1_freq)) ] )
    else:
        if j1 >= 0 and j1 < length:
            # only the first matrix
            return sum( [ c1_freq[i][j1]*math.log(c1_freq[i][j1]/bg_freq[i][0], 2) for i in range(0, len(c1_freq)) ] )
        else:
            # only the second matrix
            return sum( [ bg_freq[i][0]*math.log(bg_freq[i][0]/c2_freq[i][j2], 2) for i in range(0, len(c1_freq)) ] )

def kl_divergence(k, bg_freq, c1_freq, c2_freq):
    """ KL-divergence normalized by length """
    length  = len(c1_freq[0])
    j_range = range(min(0, k), max(length, length+k))
    return 1.0/float(length+abs(k)) * sum([ kl_divergence_column(j, k, bg_freq, c1_freq, c2_freq, length) for j in j_range ])
