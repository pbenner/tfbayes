## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## dpm library
lib_LTLIBRARIES  = libdpm.la
libdpm_la_SOURCES = data.cc data.hh dpm.cc dpm.hh clusters.cc clusters.hh statistics.cc statistics.hh init.cc init.hh interface.cc interface.hh
libdpm_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libdpm_la_LDFLAGS = -no-undefined -version-info 0:0:0

## clean python files
clean-local:
	$(RM) *.pyc
