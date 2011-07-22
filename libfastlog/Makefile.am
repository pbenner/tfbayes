## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## bunzip log table
libfastlog.c: libfastlog.c.bz2
	bunzip2 -k $<

## library
noinst_LTLIBRARIES = libfastlog.la

libfastlog_la_SOURCES = libfastlog.c
libfastlog_la_LDFLAGS = -no-undefined -version-info 0:0:0

## clean

clean-local:
	$(RM) libfastlog.c
