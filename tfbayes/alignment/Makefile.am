## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

## headers
noinst_HEADERS = alignment.hh alignment-hmm.hh sequence.hh

noinst_PROGRAMS = alignment-test

alignment_test_SOURCES = alignment-test.cc
alignment_test_LDADD   = $(top_builddir)/tfbayes/uipac/libtfbayes-uipac.la
alignment_test_LDADD  += $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la
alignment_test_LDADD  += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la

## python interface
pkgpython_LTLIBRARIES = interface.la

interface_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
interface_la_CXXFLAGS = 
interface_la_SOURCES  = \
	interface.cc
interface_la_LIBADD   = $(BOOST_PYTHON_LIB)
interface_la_LIBADD  += $(BOOST_SYSTEM_LIB)
interface_la_LIBADD  += -lm
interface_la_LIBADD  += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
interface_la_LIBADD  += $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la
interface_la_LIBADD  += $(top_builddir)/tfbayes/uipac/libtfbayes-uipac.la
interface_la_LDFLAGS  = -avoid-version
interface_la_LDFLAGS += -module

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/alignment

## compile python files
pkgpython_PYTHON = __init__.py tools.py

## clean python files
clean-local:
	$(RM) *.pyc
