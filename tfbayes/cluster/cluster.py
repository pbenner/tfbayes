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

import information

from plot        import *
from tools       import *

from ..alignment import select_subsequence
from ..dpm       import range_t, seq_index_t, dpm_subset_t

# basic class to store cluster
# ------------------------------------------------------------------------------

# A cluster is initialized with the counts of the sites that belong to
# it. The number of components gives the number of sites in the
# multiple alignment that belong to the cluster. This number does not
# consider the number of species that are present in the
# alignment. All statistics on the cluster are computed without
# considering gaps in the alignment!

# We use the following terminology:
# - counts: the raw nucleotide counts from the sequences
# - alpha : prior pseudo counts
# - motif : the posterior expectation of the nucleotides for all columns,
#           computed from counts and prior pseudocounts
# - logo  : a motif that is scaled by the discrete entropy
# - pwm   : Kullback-Leibler divergence of the motif to the background
#           distribution

class cluster_t():
    counts         = None # nucleotide counts
    counts_gap     = None # gap counts
    alpha          = None # prior pseudo counts
    alpha_gap      = None # gap prior pseudo counts
    components     = None # number of sites that belong to the cluster
    identifier     = None # some identifier (integer)
    cluster_type   = None # i.e. which baseline prior was used
    n              = None # alphabet size
    m              = None # cluster length
    def __init__(self, counts, counts_gap, alpha, alpha_gap, components, identifier, cluster_type, sites = None):
        # initialize counts, components, and identifier
        self.counts       = counts
        self.counts_gap   = counts_gap
        self.alpha        = alpha
        self.alpha_gap    = alpha_gap
        self.components   = components
        self.identifier   = identifier
        self.cluster_type = cluster_type
        self.n            = 4
        self.sites        = sites
        if not len(counts) == self.n:
            raise IOError('Counts matrix has invalid dimension.')
        self.m          = len(counts[0])
    def __getitem__(self, s):
        tr = lambda x: zip(*x)
        if not isinstance(s, slice):
            raise IOError("__getitem__() requires a slice object.")
        if s.start < 0 or s.stop-1 < s.start:
            raise IndexError("Index out of bounds")
        if s.stop - s.start > len(self.counts[0]):
            raise IndexError("Index out of bounds")
        counts     = tr(tr(self.counts)[s])
        counts_gap = self.counts_gap[s]
        alpha      = tr(tr(self.alpha )[s])
        alpha_gap  = self.alpha_gap[s]
        sites      = dpm_subset_t(self.sites.dpm_subset_tag())
        for site in self.sites:
            if not site.reverse():
                index = seq_index_t(site.index()[0], site.index()[1]+s.start)
                sites.insert(range_t(index, s.stop-s.start, site.reverse()))
            else:
                index = seq_index_t(site.index()[0], site.index()[1]-s.start)
                sites.insert(range_t(index, s.stop-s.start, site.reverse()))
        return cluster_t(counts, counts_gap, alpha, alpha_gap, self.components, self.identifier, self.cluster_type, sites = sites)
    def posterior_counts(self):
        counts = [ [ self.counts[i][j] + self.alpha[i][j] for j in range(self.m) ] for i in range(self.n) ]
        gaps   = [ self.counts_gap[i] + self.alpha_gap[i] for j in range(self.m) ]
        return counts + [gaps]
    def average_counts(self):
        return average_counts(self.counts, self.n, self.m)
    def frequencies(self):
        return frequencies(self.counts, self.n, self.m)
    def motif(self):
        return motif(self.counts, self.alpha, self.n, self.m)
    def pwm(self, bg):
        motif = self.motif()
        return pwm(motif, bg, self.n, self.m)
    def score(self, sequence):
        """Given a nucleotide sequence of the same length as the
           cluster, compute how well the pwm matches this sequence."""
        pwm = self.pwm(alpha)
        return score(pwm, sequence)
    def revcomp(self):
        counts     = revcomp(self.counts, self.n, self.m)
        counts_gap = [ self.counts_gap[-j-1] for j in range(self.m) ]
        alpha      = revcomp(self.alpha, self.n, self.m)
        alpha_gap  = [ self.alpha_gap[-j-1] for j in range(self.m) ]
        return cluster_t(counts, counts_gap, alpha, alpha_gap, self.components, self.identifier, self.cluster_type, sites = self.sites)
    def entropy(self):
        return information.entropy(self.motif())
    def r_sequence(self):
        return information.r_sequence(self.motif())
    def scan(self, sequence, threshold, skip_gaps=False):
        pwm    = self.pwm()
        length = len(pwm[0])
        for i in range(len(sequence)-length+1):
            result = select_subsequence(sequence, i, length, skip_gaps)
            if result:
                subsequence, j = result
                value = score(pwm, subsequence)
                if value > threshold:
                    yield (i,j,value)
