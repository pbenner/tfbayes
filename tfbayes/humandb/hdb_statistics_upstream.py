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

import interface

# global options
# ------------------------------------------------------------------------------

database = {
    'identifier' : None,
    'basedir'    : None,
    'db'         : None,
    'maf'        : None
    }

options = {
    'verbose'    : False,
    'pure'       : False
    }

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "hdb-statistics-seq [option]... DATABASE_CONFIG"
    print
    print "Options:"
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print


# load results from file
# ------------------------------------------------------------------------------

codons = [
    'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
    'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
];

def loadConfig(config_file):
    result = np.zeros(64)
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
        database['identifier'] = config_parser.get(section, 'identifier')
        database['database']   = config_parser.get(section, 'database')
        database['sequences']  = config_parser.get(section, 'sequences')
        fp = open(config_file+'.pkl', 'wb')
        pickle.dump(database, fp)
        fp.close()

    directory = os.path.dirname(config_file)
    interface.hdb_init('hdb-read-seq', 1)
    for seq in range(0, int(database['sequences'])):
#        print "Opening: "+str(seq+1)
        dbp = interface.hdb_open_ro(os.path.join(directory, database['database']), str(seq+1))
        result += interface.hdb_count_codons_upstream(dbp, 5000)
        interface.hdb_close(dbp)
    interface.hdb_free()

    # print statistics
    tmp   = {}
    total = float(sum(result))
    for i in range(0,64):
        tmp[codons[i]] = result[i]/total
    print tmp

# main
# ------------------------------------------------------------------------------

def main():
    global options
    try:
        longopts   = ['help', 'verbose']
        opts, tail = getopt.getopt(sys.argv[1:], "p", longopts)
    except getopt.GetoptError:
        usage()
        return 2
    output = None
    for o, a in opts:
        if o in ('-v', '--verbose'):
            sys.stderr.write('Verbose mode turned on.\n')
            options['verbose'] = True
        if o in ('-h', '--help'):
            usage()
            return 0
        if o in ('-p', '--pure'):
            options['pure'] = True
    if len(tail) != 1:
        usage()
        return 1
    loadConfig(tail[0])
    return 0

if __name__ == '__main__':
    sys.exit(main())
