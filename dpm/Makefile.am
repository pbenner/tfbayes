## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

bin_PROGRAMS = tfbs-dpm

tfbs_dpm_SOURCES = tfbs-dpm.cc
tfbs_dpm_LDADD   = $(top_builddir)/libexception/libexception.la
tfbs_dpm_LDADD  += libdpm.la
tfbs_dpm_LDFLAGS = -static

## dpm library
lib_LTLIBRARIES  = libdpm.la
libdpm_la_SOURCES = data.cc data.hh dpm.cc dpm.hh clustermanager.cc clustermanager.hh cluster.cc cluster.hh distribution.cc distribution.hh sampler.cc sampler.hh init.cc init.hh interface.cc interface.hh
libdpm_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libdpm_la_LDFLAGS = -no-undefined -version-info 0:0:0

## clean python files
clean-local:
	$(RM) *.pyc
