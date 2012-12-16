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

from ..uipac.alphabet import DNA

# counts -> frequencies
# ------------------------------------------------------------------------------

def frequencies(counts, n, m):
    sums =   [ sum(map(lambda c: float(c[j]), counts)) for j in range(m) ]
    freq = [ [ float(counts[i][j])/sums[j]             for j in range(m) ] for i in range(n) ]
    return freq

# score of a sequence given a pwm
# ------------------------------------------------------------------------------

def score(pwm, sequence):
    result = 0
    for j in range(len(sequence)):
        if sequence[j] == 'A' or sequence[j] == 'a':
            code = DNA.code('A')
            result += pwm[code][j]
        if sequence[j] == 'C' or sequence[j] == 'c':
            code = DNA.code('C')
            result += pwm[code][j]
        if sequence[j] == 'G' or sequence[j] == 'g':
            code = DNA.code('G')
            result += pwm[code][j]
        if sequence[j] == 'T' or sequence[j] == 't':
            code = DNA.code('T')
            result += pwm[code][j]
    return result

# reverse complement
# ------------------------------------------------------------------------------

def revcomp(counts, n, m):
    result = [ [ 0.0 for j in range(m) ] for i in range(n) ]
    for j in range(m):
        result[DNA.code('T')][-j-1] = counts[DNA.code('A')][j]
        result[DNA.code('G')][-j-1] = counts[DNA.code('C')][j]
        result[DNA.code('C')][-j-1] = counts[DNA.code('G')][j]
        result[DNA.code('A')][-j-1] = counts[DNA.code('T')][j]
    return result

# average number of nucleotides a logo is constructed from
# ------------------------------------------------------------------------------

def average_counts(counts, n, m):
    return sum([ sum([ counts[i][j] for i in range(n) ]) for j in range(m) ])/float(m)

# compute the posterior expectations of a cluster given
# prior pseudo counts
# ------------------------------------------------------------------------------

def motif(counts, alpha, n, m):
    counts_sum  = [ sum(map(lambda c: float(c[j]), counts)) for j in range(m) ]
    alpha_sum   = [ sum(map(lambda a: float(a[j]), alpha )) for j in range(m) ]
    expectation = [ [ float(counts[i][j] + alpha[i][j])/(counts_sum[j]+alpha_sum[j])
                      for j in range(m) ]
                    for i in range(n) ]
    return expectation

# motif -> pwm
# ------------------------------------------------------------------------------

def pwm(motif, n, m):
    bg = [ 0.0 ] * n
    bg[DNA.code('A')] = 0.3
    bg[DNA.code('C')] = 0.2
    bg[DNA.code('G')] = 0.2
    bg[DNA.code('T')] = 0.3
    return [ [ math.log(motif[i][j]/bg[i], 2) for j in range(m) ] for i in range(n) ]
