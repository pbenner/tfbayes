## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## Berkeley DB for the human genome data
lib_LTLIBRARIES  = libhumandb.la
libhumandb_la_SOURCES = humandb.c humandb.h statistics-codons.c statistics-codons.h interface.c interface.h
libhumandb_la_LIBADD  = $(LIB_DB) $(LIB_PTHREAD) $(top_builddir)/libexception/libexception.la
libhumandb_la_LDFLAGS = -no-undefined -version-info 0:0:0

## set python package name
pkgpythondir = $(pyexecdir)/humandb

## compile python files
pkgpython_PYTHON = __init__.py hdb_load_upstream.py hdb_load_chrom.py hdb_read_upstream.py hdb_read_chrom.py hdb_statistics_chrom.py hdb_statistics_upstream.py hdb_search_upstream.py hdb_search_chrom.py config.py interface.py statistics.py

## clean python files
clean-local:
	$(RM) *.pyc
