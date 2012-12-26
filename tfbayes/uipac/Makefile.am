## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)

## library
lib_LTLIBRARIES = libtfbayes-uipac.la

libtfbayes_uipac_la_SOURCES = \
	code.cc \
	code.hh \
	nucleotide-sequence.hh
libtfbayes_uipac_la_LDFLAGS = -no-undefined -version-info 0:0:0

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/uipac

## compile python files
pkgpython_PYTHON = __init__.py alphabet.py

## clean python files
clean-local:
	$(RM) *.pyc
