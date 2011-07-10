#! /usr/bin/env python

# imports
# ------------------------------------------------------------------------------

import re
import math
import time
import gobject
import ConfigParser
import gtk
import sys
import getopt
import os

import numpy               as     np
import numpy.random.mtrand as     mt
import numpy.random        as     rd

def importMatplotlib(backend=None):
    global pyplot
    global NonUniformImage
    global PdfPages
    global patches
    global path
    from matplotlib import use
    if backend:
        use(backend)
    import matplotlib.pyplot   as     pyplot
    from   matplotlib.image    import NonUniformImage
    from   matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.patches  as     patches
    import matplotlib.path     as     path

importMatplotlib('Agg')

# local imports
# ------------------------------------------------------------------------------

from interface import *

# global options
# ------------------------------------------------------------------------------

options = {
    'save'        : None,
    'interactive' : False,
    'verbose'     : False,
    }

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "tfbs-plot.py [option]... RESULT SEQUENCES"
    print
    print "Options:"
    print "   -s, --save=FILE                - save plot to FILE"
    print
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print

# parse sequences file
# ------------------------------------------------------------------------------

def makefilter(keep):
    """ Return a functor that takes a string and returns a copy of that
        string consisting of only the characters in 'keep'.
    """
    import string

    # make a string of all chars, and one of all those NOT in 'keep'
    allchars = string.maketrans('', '')
    delchars = ''.join([c for c in allchars if c not in keep])

    # return the functor
    return lambda s,a=allchars,d=delchars: s.translate(a, d)

def readSequences(seq_file):
    fh = open(seq_file,'r')
    identifier = makefilter("ACGTacgt")
    li = fh.readlines()
    sequences = []
    for line in li:
        if line != "\n":
            sequences.append(identifier(line))
    fh.close()
    return sequences

# parse results config
# ------------------------------------------------------------------------------

def readVector(config, section, option, converter):
    vector_str = config.get(section, option)
    vector     = map(converter, vector_str.split(' '))
    return vector

def readMatrix(config, section, option, converter):
    matrix_str = config.get(section, option)
    matrix     = []
    for line in matrix_str.split('\n'):
        if line != '':
            matrix.append([converter(a) for a in line.split(' ')])
    return matrix

def parseConfig(config_file):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)

    posterior  = readMatrix(config_parser, 'Result', 'posterior', float)
    likelihood = readVector(config_parser, 'Result', 'likelihood', float)
    components = readVector(config_parser, 'Result', 'components', int)

    return posterior, likelihood, components

# plot
# ------------------------------------------------------------------------------

def plot_likehood(likelihood, pp):
    fig = pyplot.figure()
    ax  = fig.add_subplot(111, title="Log Likelihood")
    ax.plot(likelihood)
    pp.savefig()

def plot_components(components, pp):
    fig = pyplot.figure()
    ax  = fig.add_subplot(111, title="Components")
    ax.plot(components)
    pp.savefig()

def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def plot_sequence(sequence, _posterior, n, pp):
    sequences = split_len(sequence, 50)
    posterior = split_len(_posterior, 50)
    fig = pyplot.figure()
    ax  = fig.add_subplot(111, title="Sequence "+str(n))
    ax.set_frame_on(False)
    ax.set_axis_off()
    xn    = max(map(len, sequences))
    yn    = len(sequences)
    yfrom = 0.90
    yto   = yfrom-yn*0.07
    x, y = np.meshgrid(np.arange(0.01, 1.00,  0.99/xn),
                       np.arange(yfrom, yto, -(yfrom-yto)/len(sequences)))
    for sl, po, p1l, p2l in zip(sequences, posterior, x, y):
        for s, p, p1, p2 in zip(sl, po, p1l, p2l):
            ax.text(p1, p2, str(s), size=10, rotation=0,
                    ha="center", va="bottom", color=pyplot.cm.jet(p))
    pp.savefig()

def plot_sequences(sequences, posterior, pp):
    for i in range(0, len(sequences)):
        plot_sequence(sequences[i], posterior[i], i, pp)

def plot_result(results_file, sequences_file):
    posterior, likelihood, components = parseConfig(results_file)
    sequences = readSequences(sequences_file)
    pp = PdfPages(options['save'])
    plot_likehood(likelihood, pp)
    plot_components(components, pp)
    plot_sequences(sequences, posterior, pp)
    pp.close()

# main
# ------------------------------------------------------------------------------

def main():
    global options
    try:
        longopts   = ["help", "verbose", "save="]
        opts, tail = getopt.getopt(sys.argv[1:], "s:", longopts)
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
        if o in ("-s", "--save"):
            options['save'] = a
    if len(tail) != 2:
        usage()
        return 1
    plot_result(tail[0], tail[1])
    return 0

if __name__ == "__main__":
    sys.exit(main())
