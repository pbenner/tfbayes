#! /usr/bin/env python

# imports
# ------------------------------------------------------------------------------

import re
import math
import time
import gobject
import gtk
import sys
import getopt
import os

import numpy               as     np
import numpy.random.mtrand as     mt
import numpy.random        as     rd

from   matplotlib          import *
use('GTKAgg')

from   matplotlib.pyplot   import *
from   matplotlib.image    import NonUniformImage
import matplotlib.patches  as     patches
import matplotlib.path     as     path

# local imports
# ------------------------------------------------------------------------------

from interface  import *

# global options
# ------------------------------------------------------------------------------

options = {
    'verbose' : False
    }

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "hdb-load-seq [option]... DATABASE_CONFIG INPUT_FILE SEQUENCE_NUMBER"
    print
    print "Options:"
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print

# Tfbs DPM
# ------------------------------------------------------------------------------

class TfbsDPM():
    def __init__(self, sequences, classes):
        dpm_init(sequences, classes)
    def print_clusters(self):
        dpm_print()
    def num_clusters(self):
        return dpm_num_clusters()
    def sample(self, n):
        dpm_sample(n)
        self.steps += 10
    def plotData(self, ax):
        return None
    def plotResult(self, ax):
        return None
    def plotStatistics(self, ax1):
        ax2 = ax1.twinx()
        p1  = ax1.plot(dpm_hist_switches())
        p2  = ax2.plot(dpm_hist_likelihood(), color='green')
        ax1.set_ylabel("Class switches")
        ax2.set_ylabel("Likelihood")
        ax1.set_xlabel("iteration")
        ax1.legend([p1, p2], ["Mean class switches", "Mean likelihood"])

class InteractiveTDPM(TfbsDPM):
    def __init__(self, sequences, classes, ax):
        TfbsDPM.__init__(self, sequences, classes)
        self.plotResult(ax)
        manager = get_current_fig_manager()
        def updatefig(*args):
            try:
                ax.cla()
                self.sampleInteractively(1, ax)
                self.plotResult(ax)
                manager.canvas.draw()
                return True
            except StopIteration:
                return False
        gobject.idle_add(updatefig)
    def sampleInteractively(self, n, ax):
#        self.sample(n)
        return None

# main
# ------------------------------------------------------------------------------

def readSequences(seq_file):
    fh = open(seq_file,'r')
    li = fh.readlines()
    sequences = []
    classes   = []
    for line in li:
        sequences.append('')
        classes.append([])
        for m in re.finditer(r'\(([1-9]):([ACGT]+)\)|[acgt]+', line):
            if m.group(2):
                classes[-1].extend([ int(m.group(1)) for i in range(0, len(m.group(2))) ])
                sequences[-1] = sequences[-1] + m.group(2)
            else:
                classes[-1].extend([ 0 for i in range(0, len(m.group(0))) ])
                sequences[-1] = sequences[-1] + m.group(0)
    fh.close()
    return sequences, classes

def sample(sequences, classes):
    fig1  = figure()
    ax1   = fig1.add_subplot(2,1,1, title="Data")
    ax2   = fig1.add_subplot(2,1,2)
    dpm   = InteractiveTDPM(sequences, classes, ax2)
    dpm.plotData(ax1)
    dpm.plotResult(ax2)
    show()

    fig2 = figure()
    ax3  = fig2.add_subplot(1,1,1, title="Statistics")
    dpm.plotStatistics(ax3)
    show()

def main():
    global options
    try:
        longopts   = ["help", "verbose"]
        opts, tail = getopt.getopt(sys.argv[1:], "", longopts)
    except getopt.GetoptError:
        usage()
        return 2
    output = None
    for o, a in opts:
        if o in ("-v", "--verbose"):
            sys.stderr.write("Verbose mode turned on.\n")
            options["verbose"] = True
        if o in ("-h", "--help"):
            usage()
            return 0
    if len(tail) != 1:
        usage()
        return 1
    sequences, classes = readSequences(tail[0])
    sample(sequences, classes)
    return 0

if __name__ == "__main__":
    sys.exit(main())
