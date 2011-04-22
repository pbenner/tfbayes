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

verbose = False

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "hdb-search-seq [option]... DATABASE_BUNDLE SEQUENCE SEQUENCE_NUMBER"
    print
    print "Options:"
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print


# load results from file
# ------------------------------------------------------------------------------

def loadConfig(config_file, sequence, seq_num):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)
    if not config_parser.has_section('Database Bundle'):
        raise IOError("Invalid configuration file.")

    databases = config.readMatrix(config_parser, 'Database Bundle', 'databases', str)
    dbp_list  = []

    interface.hdb_init('hdb-search-seq')
    for database_id, database_file in databases:
        dbp_list.append(interface.hdb_open(database_file, seq_num))

    interface.hdb_search(dbp_list, sequence)

    for dbp in dbp_list:
        interface.hdb_close(dbp)

    interface.hdb_free()

# main
# ------------------------------------------------------------------------------

options = {
    'verbose'  : False
    }

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
    if len(tail) != 3:
        usage()
        return 1
    loadConfig(tail[0], tail[1], tail[2])
    return 0

if __name__ == "__main__":
    sys.exit(main())
