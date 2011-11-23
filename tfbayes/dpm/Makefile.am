## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

bin_PROGRAMS = dpm-tfbs-debug dpm-gaussian-debug

dpm_gaussian_debug_SOURCES = dpm-gaussian-main.cc
dpm_gaussian_debug_LDADD   = $(top_builddir)/libexception/libexception.la
dpm_gaussian_debug_LDADD  += $(top_builddir)/parsmm/libparsmm.la
dpm_gaussian_debug_LDADD  += $(top_builddir)/tfbayes/fasta/libfasta.la
dpm_gaussian_debug_LDADD  += libdpm.la
dpm_gaussian_debug_LDADD  += $(LIB_PTHREAD)

dpm_tfbs_debug_SOURCES = dpm-tfbs-main.cc
dpm_tfbs_debug_LDADD   = $(top_builddir)/libexception/libexception.la
dpm_tfbs_debug_LDADD  += $(top_builddir)/parsmm/libparsmm.la
dpm_tfbs_debug_LDADD  += $(top_builddir)/tfbayes/fasta/libfasta.la
dpm_tfbs_debug_LDADD  += libdpm.la
dpm_tfbs_debug_LDADD  += $(LIB_PTHREAD)

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/dpm

## compile python files
pkgpython_PYTHON = __init__.py dpm_gaussian_interface.py dpm_tfbs_interface.py

## scripts
dist_bin_SCRIPTS = dpm-gaussian-demo dpm-tfbs-align dpm-tfbs-demo dpm-tfbs-graph dpm-tfbs-plot dpm-tfbs-sampler dpm-tfbs-scan dpm-tfbs-tools dpm-tfbs-weblogo

## dpm library
lib_LTLIBRARIES  = libdpm.la
libdpm_la_SOURCES =                                                           \
	code.cc code.hh                                                       \
	index.hh indexer.hh clonable.hh observer.hh statistics.hh             \
	datatypes.cc datatypes.hh graph.hh                                    \
	data.hh data-tfbs.cc data-tfbs.hh data-gaussian.cc data-gaussian.hh   \
	dpm-tfbs.cc dpm-tfbs.hh dpm-gaussian.cc dpm-gaussian.hh               \
	clustermanager.cc clustermanager.hh                                   \
	cluster.cc cluster.hh                                                 \
	component-model.cc component-model.hh                                 \
	mixture-model.hh                                                      \
	pmcmc.cc pmcmc.hh                                                     \
	sampler.cc sampler.hh                                                 \
	init.cc init.hh                                                       \
	interface.cc interface.hh                                             \
	dpm-gaussian-interface.cc dpm-gaussian-interface.hh                   \
	dpm-tfbs-interface.cc dpm-tfbs-interface.hh
libdpm_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libdpm_la_LIBADD += $(top_builddir)/parsmm/libparsmm.la
libdpm_la_LIBADD += $(top_builddir)/tfbayes/fasta/libfasta.la
libdpm_la_LDFLAGS = -no-undefined -version-info 0:0:0

## extra dist
EXTRA_DIST = dpm-tfbs-demo-1.cfg dpm-tfbs-demo-1.fa dpm-tfbs-demo-2.cfg dpm-tfbs-demo-2.fa

## clean python files
clean-local:
	$(RM) *.pyc
