## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

BUILT_SOURCES = libfastlog.c

## bunzip log table
libfastlog.c: libfastlog.c.bz2
	bunzip2 -k $<

## library
noinst_LTLIBRARIES = libfastlog.la

nodist_libfastlog_la_SOURCES = libfastlog.c
libfastlog_la_LDFLAGS = -no-undefined -version-info 0:0:0

## extra dist
EXTRA_DIST = libfastlog.c.bz2

## clean
CLEANFILES = libfastlog.c
