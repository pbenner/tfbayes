## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## Berkeley DB for the human genome data
lib_LTLIBRARIES  = libfasta.la
libfasta_la_SOURCES = fasta.cc fasta.hh
libfasta_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libfasta_la_LDFLAGS = -no-undefined -version-info 0:0:0

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/fasta

## compile python files
pkgpython_PYTHON = __init__.py fasta.py

## clean python files
clean-local:
	$(RM) *.pyc
