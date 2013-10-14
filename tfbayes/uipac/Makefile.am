## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/include
#AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

noinst_PROGRAMS = alphabet-test

alphabet_test_SOURCES = alphabet-test.cc
alphabet_test_LDADD   = libtfbayes-uipac.la

## library
lib_LTLIBRARIES = libtfbayes-uipac.la

libtfbayes_uipac_la_SOURCES = \
	alphabet.cc
	alphabet.hh
	code.cc \
	code.hh \
	nucleotide-sequence.hh
libtfbayes_uipac_la_LDFLAGS = -no-undefined -version-info 0:0:0

otherincludedir = $(includedir)/tfbayes/uipac
otherinclude_HEADERS = \
	code.hh \
	nucleotide-sequence.hh

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/uipac

## compile python files
pkgpython_PYTHON = __init__.py alphabet.py

## clean python files
clean-local:
	$(RM) *.pyc
