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
    'verbose'    : False,
    'pure'       : False
    }

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "hdb-read-seq [option]... DATABASE_CONFIG SEQUENCE_NUMBER POSITION NUCLEOTIDES"
    print
    print "Options:"
    print "   -p, --pure                     - strip all alignment specific characters"
    print
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print


# load results from file
# ------------------------------------------------------------------------------

def loadConfig(config_file, seq, pos, num):
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

    if not int(seq) > 0 and int(seq) > database['sequences']:
        raise IOError("Invalid sequence number.")

    interface.hdb_init('hdb-read-seq')
    dbp = interface.hdb_open(database['database'], seq)
    if options['pure']:
        print interface.hdb_get_sequence_pure(dbp, pos, num)
    else:
        print interface.hdb_get_sequence(dbp, pos, num)
    interface.hdb_close(dbp)
    interface.hdb_free()

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
    if len(tail) != 4:
        usage()
        return 1
    loadConfig(tail[0], tail[1], int(tail[2]), int(tail[3]))
    return 0

if __name__ == '__main__':
    sys.exit(main())
