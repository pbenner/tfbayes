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

import sys
import getopt
import os
import ConfigParser
import numpy as np
import math
import random
import pickle

import humandb.config    as config
import humandb.interface as interface

# global options
# ------------------------------------------------------------------------------

database = {
    'identifier' : None,
    'basedir'    : None,
    'db'         : None,
    'maf'        : None
    }

options = {
    'config-name'  : None,
    'cluster-name' : None,
    'sequence'     : None,
    'threshold'    : 20,
    'verbose'      : False
    }

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "hdb-search-chrom [option]... DATABASE_CONFIG"
    print
    print "Options:"
    print "   -p=CLUSTER_NAME:CLUSTER.cfg    - search for pwm match"
    print "   -s=SEQUENCE                    - search for exact match of a sequence"
    print
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print

# compute pwm
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

def compute_frequencies(counts):
    sums = [ sum(map(lambda m: m[j], counts)) for j in range(0, len(counts[0])) ]
    return [ [ float(counts[i][j])/sums[j]    for j in range(0, len(counts[0])) ] for i in range(0, len(counts)) ]

def compute_pwm(config_parser, cluster_name):
    bg_counts   = readMatrix(config_parser, 'Cluster', 'cluster_0',  int)
    tfbs_counts = readMatrix(config_parser, 'Cluster', cluster_name, int)
    bg_freq     = compute_frequencies(bg_counts)
    tfbs_freq   = compute_frequencies(tfbs_counts)
    return [ [ math.log(tfbs_freq[i][j]/bg_freq[i][0], 2) for j in range(0, len(tfbs_freq[0])) ] for i in range(0, len(tfbs_freq)) ]

# search
# ------------------------------------------------------------------------------

def search_pwm(dbp_list, database):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(options['config-name'])
    if not config_parser.has_section('Cluster'):
        raise IOError("Invalid configuration file.")
    pwm = compute_pwm(config_parser, options['cluster-name'])
    interface.hdb_search_pwm(dbp_list, database['chromosomes'], pwm, options['threshold'])

def search_sequence(dbp_list, database):
    interface.hdb_search(dbp_list, database['chromosomes'], options['sequence'])

# load results from file
# ------------------------------------------------------------------------------

def loadConfig(config_file):
    global database
    if os.path.isfile(config_file+'.pkl'):
        fp = open(config_file+'.pkl', 'r')
        database = pickle.load(fp)
        fp.close()
    else:
        config_parser = ConfigParser.RawConfigParser()
        config_parser.read(config_file)
        if not config_parser.has_section('Database'):
            raise IOError("Invalid configuration file.")
        section = 'Database'
        database['identifier']   = config_parser.get(section, 'identifier')
        database['database']     = config_parser.get(section, 'database')
        database['chromosomes']  = config_parser.get(section, 'chromosomes').split(' ')
        fp = open(config_file+'.pkl', 'wb')
        pickle.dump(database, fp)
        fp.close()

    directory = os.path.dirname(config_file)
    interface.hdb_init('hdb-search-chrom')
    dbp_list = []
    for chrom_name in database['chromosomes']:
        dbp_list.append(interface.hdb_open_ro(os.path.join(directory, database['database']), chrom_name))

    if options['sequence']:
        search_sequence(dbp_list, database)
    if options['cluster-name']:
        search_pwm(dbp_list, database)

    for dbp in dbp_list:
        interface.hdb_close(dbp)

    interface.hdb_free()

# main
# ------------------------------------------------------------------------------

def main():
    global options
    try:
        longopts   = ["help", "verbose", "threshold="]
        opts, tail = getopt.getopt(sys.argv[1:], "p:s:", longopts)
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
        if o == "-p":
            options['config-name'], options['cluster-name'] = a.split(':')
            if options['config-name'] == '' or options['cluster-name'] == '':
                usage()
                return 1
        if o == "-s":
            options['sequence'] = a
        if o == "--threshold":
            options['threshold'] = float(a)
    if len(tail) != 1:
        usage()
        return 1
    loadConfig(*tail)
    return 0

if __name__ == "__main__":
    sys.exit(main())
