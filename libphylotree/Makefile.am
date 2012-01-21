## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

bin_PROGRAMS = polynomial-test

polynomial_test_SOURCES = polynomial-test.cc
polynomial_test_LDADD   = -lm

noinst_LTLIBRARIES = libphylotree.la

libphylotree_la_SOURCES = phylotree.cc phylotree.hh
libphylotree_la_LIBADD  = -lm
libphylotree_la_LDFLAGS = -no-undefined -version-info 0:0:0
