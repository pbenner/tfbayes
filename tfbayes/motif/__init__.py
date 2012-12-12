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

from plot import *

# tools
# ------------------------------------------------------------------------------

def average_components(motif):
    counts  = [ sum([ motif[j][i] for j in range(4) ]) for i in range(len(motif[0])) ]
    average = sum(counts)/float(len(motif[0]))
    return average

# 
# ------------------------------------------------------------------------------

class cluster_t():
    counts         = None
    average_counts = None
    components     = None
    identifier     = None
    def __init__(self, n, m):
        # initialize counts with a simple matrix of zeros
        self.counts = [ [ 0.0 for j in range(n) ] for i in range(m) ]
