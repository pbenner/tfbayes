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
import re

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
    'chromosome' : None,
    'gene_type'  : None,
    'verbose'    : False,
    'pure'       : False
    }

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "hdb-read-chrom [option]... DATABASE_CONFIG ENCODE_GTF_FILE"
    print
    print "Options:"
    print "   --chromosome=CHROMOSOME        - use only the specified chromosome"
    print "   --gene-type=GENE_TYPE          - use only the specified gene type"
    print
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print

# parse encode file
# ------------------------------------------------------------------------------

codons = [
    'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
    'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
];

def parseDescription(description):
    tmp = re.split('; |;', description)
    description = {}
    for entry in tmp:
        if not entry == '' and not entry == '\n':
            key, value = re.split(' ', entry, maxsplit=1)
            description[key] = value.replace('\"','')
    return description

def parseFile(filename, directory, database):
    fp = open(filename, 'r')
    positions = {}
    result = np.zeros(64)

    # this parses .gtf files from the encode project
    for line in fp.readlines():
        s = re.split('\t', line)
        if (len(s) == 9 and (not options['chromosome'] or s[0] == options['chromosome'])):
            description = parseDescription(s[8])
            if not options['gene_type'] or description['gene_type'] == options['gene_type']:
                try:
                    positions[s[0]]
                except KeyError:
                    positions[s[0]] = []
                positions[s[0]].append(long(s[3]))
                positions[s[0]].append(long(s[4]))
    fp.close()

    # for each chromosome, count the number of codons at the positions
    # specified in the array 'positions'
    for chrom in positions:
        try:
            database['chromosomes'].index(chrom)
        except ValueError:
            print "Unknown chromosome."
            exit(1)
        dbp = interface.hdb_open_ro(os.path.join(directory, database['database']), chrom)
        result += interface.hdb_count_codons(dbp, positions[chrom])
        interface.hdb_close(dbp)

    # print statistics
    tmp   = {}
    total = float(sum(result))
    for i in range(0,64):
        tmp[codons[i]] = result[i]/total
    print tmp

# load results from file
# ------------------------------------------------------------------------------

def loadConfig(config_file, filename):
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
    interface.hdb_init('hdb-statistics-chrom', 1)

    parseFile(filename, directory, database)

    interface.hdb_free()

# main
# ------------------------------------------------------------------------------

def main():
    global options
    try:
        longopts   = ['help', 'verbose', 'chromosome=', 'gene-type=']
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
        if o == '--chromosome':
            options['chromosome'] = a
        if o == '--gene-type':
            options['gene_type'] = a
        if o in ('-p', '--pure'):
            options['pure'] = True
    if len(tail) == 2:
        loadConfig(tail[0], tail[1])
    else:
        usage()
        return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
