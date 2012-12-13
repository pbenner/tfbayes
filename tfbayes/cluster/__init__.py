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

import fileio
import plot
import tools

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
    n              = None # alphabet size
    m              = None # cluster length
    def __init__(self, counts, counts_gap, alpha, alpha_gap, components, identifier):
        # initialize counts, components, and identifier
        self.counts     = counts
        self.counts_gap = counts_gap
        self.alpha      = alpha
        self.alpha_gap  = alpha_gap
        self.components = components
        self.identifier = identifier
        self.n          = 4
        if not len(counts) == self.n:
            raise ValueError('Counts matrix has invalid dimension.')
        self.m          = len(counts[0])
    def average_counts(self):
        return tools.average_counts(self.counts, self.n, self.m)
    def frequencies(self):
        return tools.frequencies(self.counts, self.n, self.m)
    def motif(self):
        return tools.motif(self.counts, self.alpha, self.n, self.m)
    def pwm(self):
        motif = self.motif()
        return tools.pwm(motif, self.n, self.m)
    def score(self, sequence):
        """Given a nucleotide sequence of the same length as the
           cluster, compute how well the pwm matches this sequence."""
        pwm = self.pwm(alpha)
        return tools.score(pwm, sequence)
    def revcomp(self):
        counts     = tools.revcomp(self.counts, self.n, self.m)
        counts_gap = [ self.counts_gap[-j-1] for j in range(self.m) ]
        alpha      = tools.revcomp(self.alpha, self.n, self.m)
        alpha_gap  = [ self.alpha_gap[-j-1] for j in range(self.m) ]
        return cluster_t(counts, counts_gap, alpha, alpha_gap, self.components, self.identifier)
    def entropy(self):
        motif = self.motif()
        return information.entropy(motif, self.n, self.m)
    def r_sequence(self):
        motif = self.motif()
        return information.r_sequence(motif, self.n, self.m)
