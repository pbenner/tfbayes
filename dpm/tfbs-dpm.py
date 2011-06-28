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

from interface import *

# global options
# ------------------------------------------------------------------------------

options = {
    'samples'     : (100,1000),
    'save'        : None,
    'interactive' : False,
    'verbose'     : False,
    }

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "tfbs-dpm.py [option]... SEQUENCES"
    print
    print "Options:"
    print "   -i                             - interactive mode"
    print "   -s, --save=FILE                - save posterior to FILE"
    print "       --samples=BURN_IN:SAMPLES  - number of samples [default: 100:1000]"
    print
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print

# Tfbs DPM
# ------------------------------------------------------------------------------

def cluster_colors(c):
    if c < 5:
        return {
           -1 : 'black',
            0 : 'black',
            1 : 'red',
            2 : 'green',
            3 : 'blue',
            4 : 'yellow'
            }[c]
    else:
        return 'magenta'

def plot_sequences(ax, sequences, clusters):
    ax.set_frame_on(False)
    xn    = max(map(len, sequences))
    yn    = len(sequences)
    yfrom = 0.90
    yto   = yfrom-yn*0.07
    x, y = np.meshgrid(np.arange(0.01, 1.00,  0.99/xn),
                       np.arange(yfrom, yto, -(yfrom-yto)/len(sequences)))
    for sl, cl, p1l, p2l in zip(sequences, clusters, x, y):
        for s, c, p1, p2 in zip(sl, cl, p1l, p2l):
            ax.text(p1, p2, str(s), size=10, rotation=0,
                    ha="center", va="bottom", color=cluster_colors(c))

class TfbsDPM():
    def __init__(self, sequences, clusters):
        self.sequences = sequences
        self.clusters  = clusters
        self.steps     = 0
        dpm_init(sequences[:][:], clusters[:][:])
    def print_clusters(self):
        dpm_print()
    def num_clusters(self):
        return dpm_num_clusters()
    def sample(self, n, burnin):
        dpm_sample(n, burnin)
        self.steps += n
    def plotData(self, ax):
        ax.set_xticks([]); ax.set_yticks([])
        plot_sequences(ax, self.sequences, self.clusters)
    def plotResult(self, ax):
        ax.set_title('Gibbs Sampling')
        ax.set_xticks([]); ax.set_yticks([])
        sequences = np.array(self.sequences)
        clusters = []
        for i in self.clusters:
            clusters.append(-np.ones_like(i))
        num_clusters = dpm_num_clusters()
        for c in range(0, num_clusters):
            for seq, pos in dpm_cluster(c):
                clusters[int(seq)][int(pos)] = c
        plot_sequences(ax, sequences, clusters)
    def plotPosterior(self, ax):
        ax.set_title('Posterior')
        ax.set_xticks([]); ax.set_yticks([])
        X = dpm_get_posterior()
        ax.imshow(X, cmap=cm.jet, interpolation='nearest')
    def plotStatistics(self, ax1):
        ax2 = ax1.twinx()
        p1  = ax1.plot(dpm_hist_switches())
        p2  = ax2.plot(dpm_hist_likelihood(), color='green')
        ax1.set_ylabel("Class switches")
        ax2.set_ylabel("Likelihood")
        ax1.set_xlabel("iteration")
        ax1.legend([p1, p2], ["Mean class switches", "Mean likelihood"])

class InteractiveTDPM(TfbsDPM):
    def __init__(self, sequences, clusters, ax2, ax3):
        TfbsDPM.__init__(self, sequences, clusters[:][:])
        self.plotResult(ax2)
        self.plotPosterior(ax3)
        manager = get_current_fig_manager()
        def updatefig(*args):
            try:
                ax2.cla()
                ax3.cla()
                self.sampleInteractively(1)
                self.plotResult(ax2)
                self.plotPosterior(ax3)
                manager.canvas.draw()
                return True
            except StopIteration:
                return False
        gobject.idle_add(updatefig)
    def sampleInteractively(self, n):
        self.sample(n, 0)

# main
# ------------------------------------------------------------------------------

def readSequences(seq_file):
    fh = open(seq_file,'r')
    li = fh.readlines()
    sequences = []
    clusters  = []
    for line in li:
        sequences.append('')
        clusters.append([])
        for m in re.finditer(r'\(([1-9]):([ACGT]+)\)|[ACGT]+', line):
            if m.group(2):
                clusters[-1].extend([ int(m.group(1)) for i in range(0, len(m.group(2))) ])
                sequences[-1] = sequences[-1] + m.group(2)
            else:
                clusters[-1].extend([ 0 for i in range(0, len(m.group(0))) ])
                sequences[-1] = sequences[-1] + m.group(0)
    fh.close()
    return sequences, clusters

def save_posterior(dpm):
    fp = open(options['save'], 'w')
    fp.write(str(dpm_get_posterior()))
    fp.close()

def sample(sequences, clusters):
    fig1  = figure(linewidth=0,facecolor='white',edgecolor='white')
    ax1   = fig1.add_subplot(3,1,1, title='Sequences')
    ax2   = fig1.add_subplot(3,1,2, title='Gibbs Sampling')
    ax3   = fig1.add_subplot(3,1,3, title='Posterior')
    if options['interactive']:
        dpm = InteractiveTDPM(sequences, clusters[:][:], ax2, ax3)
        dpm.plotData(ax1)
        dpm.sample(0, options['samples'][0])
        dpm.sampleInteractively(options['samples'][1])
        show()
    else:
        dpm = TfbsDPM(sequences, clusters[:][:])
        dpm.sample(options['samples'][1], options['samples'][0])
        if options['save']:
            save_posterior(dpm)
        else:
            dpm.plotData(ax1)
            dpm.plotResult(ax2)
            dpm.plotPosterior(ax3)
            show()

def main():
    global options
    try:
        longopts   = ["help", "verbose", "save=", "samples="]
        opts, tail = getopt.getopt(sys.argv[1:], "is:", longopts)
    except getopt.GetoptError:
        usage()
        return 2
    output = None
    for o, a in opts:
        if o in ("-v", "--verbose"):
            sys.stderr.write("Verbose mode turned on.\n")
            options["verbose"] = True
        if o == "-i":
            options["interactive"] = True
        if o in ("-h", "--help"):
            usage()
            return 0
        if o in ("-s", "--save"):
            options['save'] = a
        if o == "--samples":
            options["samples"] = tuple(map(int, a.split(":")))
    if len(tail) != 1:
        usage()
        return 1
    sequences, clusters = readSequences(tail[0])
    sample(sequences, clusters)
    return 0

if __name__ == "__main__":
    sys.exit(main())
