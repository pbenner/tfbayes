## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)

## install headers
otherincludedir = $(includedir)/tfbayes/alignment
otherinclude_HEADERS = alignment.hh alignment-hmm.hh

noinst_PROGRAMS = alignment-hmm-test

alignment_hmm_test_SOURCES = alignment-hmm-test.cc
alignment_hmm_test_LDADD   = -lm
alignment_hmm_test_LDADD  += $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la
alignment_hmm_test_LDADD  += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la

## library
lib_LTLIBRARIES = libtfbayes-alignment.la

libtfbayes_alignment_la_SOURCES = \
	interface.cc
libtfbayes_alignment_la_LIBADD  = -lm
libtfbayes_alignment_la_LIBADD += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
libtfbayes_alignment_la_LIBADD += $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la
libtfbayes_alignment_la_LDFLAGS = -no-undefined -version-info 0:0:0

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/alignment

## compile python files
pkgpython_PYTHON = __init__.py interface.py tools.py

## clean python files
clean-local:
	$(RM) *.pyc
