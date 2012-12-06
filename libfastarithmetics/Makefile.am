## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

bin_PROGRAMS = fast-test

fast_test_SOURCES = fast-test.cc
fast_test_LDADD   = -lm
fast_test_LDADD  += libfastarithmetics.la

## library
noinst_LTLIBRARIES = libfastarithmetics.la

libfastarithmetics_la_SOURCES  = fast-lngamma.cc
libfastarithmetics_la_SOURCES += fast-lngamma.hh
libfastarithmetics_la_SOURCES += fast-lngamma-table.cc
libfastarithmetics_la_SOURCES += fast-lngamma-table.hh
libfastarithmetics_la_SOURCES += fast-lnbeta.hh
libfastarithmetics_la_LIBADD   = -lgsl
