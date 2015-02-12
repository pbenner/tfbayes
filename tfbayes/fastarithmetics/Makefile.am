## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

noinst_PROGRAMS = fast-test

fast_test_SOURCES = fast-test.cc
fast_test_LDADD   = -lm
fast_test_LDADD  += libtfbayes-fastarithmetics.la

## library
noinst_LTLIBRARIES = libtfbayes-fastarithmetics.la

libtfbayes_fastarithmetics_la_SOURCES  = fast-lngamma.cc
libtfbayes_fastarithmetics_la_SOURCES += fast-lngamma.hh
libtfbayes_fastarithmetics_la_SOURCES += fast-lngamma-table.cc
libtfbayes_fastarithmetics_la_SOURCES += fast-lngamma-table.hh
libtfbayes_fastarithmetics_la_SOURCES += fast-lnbeta.hh
