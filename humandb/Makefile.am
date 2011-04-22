## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## Berkeley DB for the human genome data
lib_LTLIBRARIES  = libhumandb.la
libhumandb_la_SOURCES = humandb.c humandb.h interface.c interface.h
libhumandb_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libhumandb_la_LDFLAGS = -no-undefined -version-info 0:0:0

## set python package name
pkgpythondir = $(pyexecdir)/humandb

## compile python files
pkgpython_PYTHON = __init__.py hdb_load_seq.py hdb_read_seq.py hdb_search_seq.py config.py interface.py

## clean python files
clean-local:
	$(RM) *.pyc
