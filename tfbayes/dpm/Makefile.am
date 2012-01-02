## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

bin_PROGRAMS = dpm-tfbs-debug dpm-gaussian-debug parsmm-test mcm-test

parsmm_test_SOURCES = parsmm-test.cc
parsmm_test_LDADD   = $(top_builddir)/parsmm/libparsmm.la
parsmm_test_LDADD  += libdpm.la
parsmm_test_LDADD  += -lm $(LIB_PTHREAD) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

mcm_test_SOURCES = mcm-test.cc
mcm_test_LDADD   = libdpm.la
mcm_test_LDADD  += -lm $(LIB_PTHREAD) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

dpm_gaussian_debug_SOURCES = dpm-gaussian-main.cc
dpm_gaussian_debug_LDADD   = $(top_builddir)/libexception/libexception.la
dpm_gaussian_debug_LDADD  += $(top_builddir)/parsmm/libparsmm.la
dpm_gaussian_debug_LDADD  += $(top_builddir)/tfbayes/fasta/libfasta.la
dpm_gaussian_debug_LDADD  += libdpm.la
dpm_gaussian_debug_LDADD  += $(LIB_PTHREAD) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

dpm_tfbs_debug_SOURCES = dpm-tfbs-main.cc
dpm_tfbs_debug_LDADD   = $(top_builddir)/libexception/libexception.la
dpm_tfbs_debug_LDADD  += $(top_builddir)/parsmm/libparsmm.la
dpm_tfbs_debug_LDADD  += $(top_builddir)/tfbayes/fasta/libfasta.la
dpm_tfbs_debug_LDADD  += libdpm.la
dpm_tfbs_debug_LDADD  += $(LIB_PTHREAD) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/dpm

## compile python files
pkgpython_PYTHON = __init__.py dpm_gaussian_interface.py dpm_tfbs_interface.py

## dpm library
lib_LTLIBRARIES  = libdpm.la
libdpm_la_SOURCES =                                                           \
	abysmal-stack.hh code.cc code.hh                                      \
	index.cc index.hh indexer.hh clonable.hh observer.hh statistics.hh    \
	datatypes.hh graph.hh                                                 \
	data.hh data-tfbs.cc data-tfbs.hh data-gaussian.cc data-gaussian.hh   \
	dpm-tfbs.cc dpm-tfbs.hh dpm-tfbs-prior.cc dpm-tfbs-prior.hh           \
	dpm-tfbs-command.cc dpm-tfbs-command.hh                               \
	dpm-tfbs-options.hh mixture-weights.hh                                \
	dpm-tfbs-repl.cc dpm-tfbs-repl.hh                                     \
	dpm-tfbs-state.cc dpm-tfbs-state.hh dpm-tfbs-test.cc                  \
	dpm-tfbs-sampler.cc dpm-tfbs-sampler.hh                               \
	dpm-gaussian.cc dpm-gaussian.hh                                       \
	state.cc state.hh mixture-state.cc mixture-state.hh                   \
	cluster.cc cluster.hh                                                 \
	component-model.cc component-model.hh                                 \
	mixture-model.hh                                                      \
	pmcmc.cc pmcmc.hh                                                     \
	sampler.cc sampler.hh                                                 \
	init.cc init.hh utility.hh                                            \
	interface.cc interface.hh                                             \
	dpm-gaussian-interface.cc dpm-gaussian-interface.hh                   \
	dpm-tfbs-interface.cc dpm-tfbs-interface.hh
libdpm_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libdpm_la_LIBADD += $(top_builddir)/parsmm/libparsmm.la
libdpm_la_LIBADD += $(top_builddir)/tfbayes/fasta/libfasta.la
libdpm_la_LIBADD += $(LIB_PTHREAD) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)
libdpm_la_LDFLAGS = -no-undefined -version-info 0:0:0

## clean python files
clean-local:
	$(RM) *.pyc
