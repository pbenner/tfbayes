## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

## Berkeley DB for the human genome data
lib_LTLIBRARIES  = libtfbayes-fasta.la
libtfbayes_fasta_la_SOURCES = fasta.cc fasta.hh
libtfbayes_fasta_la_LIBADD  = ../exception/libtfbayes-exception.la
libtfbayes_fasta_la_LDFLAGS = -no-undefined -version-info 0:0:0

otherincludedir = $(includedir)/tfbayes/fasta
otherinclude_HEADERS = \
	fasta.hh

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/fasta

## compile python files
pkgpython_PYTHON = __init__.py fasta.py sequences.py

## clean python files
clean-local:
	$(RM) *.pyc
