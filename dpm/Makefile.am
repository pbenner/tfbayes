## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

bin_PROGRAMS = dpm-tfbs

dpm_tfbs_SOURCES = dpm-tfbs-main.cc
dpm_tfbs_LDADD   = $(top_builddir)/libexception/libexception.la
dpm_tfbs_LDADD  += libdpm.la
dpm_tfbs_LDADD  += $(LIB_PTHREAD)
dpm_tfbs_LDFLAGS = -static

## scripts
dist_bin_SCRIPTS = dpm-tfbs-plot dpm-tfbs-demo

## dpm library
lib_LTLIBRARIES  = libdpm.la
#libdpm_la_SOURCES = code.cc code.hh datatypes.cc datatypes.hh data.cc data.hh dpm.hh dpm-tfbs.cc dpm-tfbs.hh dpm-gaussian.cc dpm-gaussian.hh clustermanager.cc clustermanager.hh cluster.cc cluster.hh distribution.cc distribution.hh pmcmc.cc pmcmc.hh sampler.cc sampler.hh init.cc init.hh interface.cc interface.hh
libdpm_la_SOURCES = code.cc datatypes.cc datatypes.hh code.hh data.cc data.hh dpm.hh dpm-tfbs.cc dpm-tfbs.hh clustermanager.cc clustermanager.hh cluster.cc cluster.hh distribution.cc distribution.hh pmcmc.cc pmcmc.hh sampler.cc sampler.hh init.cc init.hh interface.cc interface.hh
libdpm_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libdpm_la_LDFLAGS = -no-undefined -version-info 0:0:0

## clean python files
clean-local:
	$(RM) *.pyc
