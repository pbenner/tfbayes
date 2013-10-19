## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/include
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

noinst_PROGRAMS = dpm-tfbs-debug dpm-gaussian-debug

# mcm_test_SOURCES = mcm-test.cc
# mcm_test_LDADD   = libtfbayes-dpm.la
# mcm_test_LDADD  += -lm $(LIB_PTHREAD) $(LIB_BOOST_REGEX) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

dpm_gaussian_debug_SOURCES = dpm-gaussian-main.cc
dpm_gaussian_debug_LDADD   = $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
dpm_gaussian_debug_LDADD  += libtfbayes-dpm.la
dpm_gaussian_debug_LDADD  += $(LIB_PTHREAD) $(LIB_BOOST_REGEX) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

dpm_tfbs_debug_SOURCES = dpm-tfbs-main.cc
dpm_tfbs_debug_LDADD   = $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
dpm_tfbs_debug_LDADD  += libtfbayes-dpm.la
dpm_tfbs_debug_LDADD  += $(LIB_PTHREAD) $(LIB_BOOST_REGEX) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

## dpm library
lib_LTLIBRARIES  = libtfbayes-dpm.la
libtfbayes_dpm_la_SOURCES =		\
	clonable.hh			\
	cluster.cc			\
	cluster.hh			\
	component-model.cc		\
	component-model.hh		\
	datatypes.hh			\
	data.hh				\
	data-tfbs.cc			\
	data-tfbs.hh			\
	data-gaussian.cc		\
	data-gaussian.hh		\
	dpm-gaussian.cc			\
	dpm-gaussian.hh			\
	dpm-sampling-history.hh		\
	dpm-tfbs.cc			\
	dpm-tfbs.hh			\
	dpm-tfbs-command.cc		\
	dpm-tfbs-command.hh		\
	dpm-tfbs-mean.cc                \
	dpm-tfbs-mean.hh                \
	dpm-tfbs-repl.cc		\
	dpm-tfbs-repl.hh		\
	dpm-tfbs-options.hh		\
	dpm-tfbs-options.cc		\
	dpm-tfbs-prior.cc		\
	dpm-tfbs-prior.hh		\
	dpm-tfbs-sampler.cc		\
	dpm-tfbs-sampler.hh		\
	dpm-tfbs-state.cc		\
	dpm-tfbs-state.hh		\
	dpm-tfbs-test.cc		\
	dpm-partition.cc		\
	dpm-partition.hh		\
	index.cc			\
	index.hh			\
	indexer.hh			\
	mixture-model.hh		\
	mixture-state.cc		\
	mixture-state.hh		\
	mixture-weights.hh		\
	observer.hh			\
	pmcmc.cc			\
	pmcmc.hh			\
	sampler.cc			\
	sampler.hh			\
	save-queue.hh                   \
	state.hh
libtfbayes_dpm_la_LIBADD  = $(top_builddir)/tfbayes/fastarithmetics/libtfbayes-fastarithmetics.la
libtfbayes_dpm_la_LIBADD += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
libtfbayes_dpm_la_LIBADD += $(top_builddir)/tfbayes/uipac/libtfbayes-uipac.la
libtfbayes_dpm_la_LIBADD += $(LIB_PTHREAD) $(LIB_BOOST_REGEX) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)
libtfbayes_dpm_la_LDFLAGS = -no-undefined -version-info 0:0:0

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/dpm

## compile python files
pkgpython_PYTHON = __init__.py

## scripts
dist_bin_SCRIPTS = \
	dpm-gaussian-demo \
	dpm-tfbs-demo

## python libraries
pkgpython_LTLIBRARIES = interface.la dpm_gaussian_interface.la dpm_tfbs_interface.la

interface_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
interface_la_CXXFLAGS =
interface_la_SOURCES  = \
	interface.cc
interface_la_LIBADD   = libtfbayes-dpm.la
interface_la_LIBADD  += $(LIB_BOOST_PYTHON)
interface_la_LDFLAGS  = -avoid-version -module

dpm_gaussian_interface_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
dpm_gaussian_interface_la_CXXFLAGS =
dpm_gaussian_interface_la_SOURCES  = \
	dpm-gaussian-interface.cc
dpm_gaussian_interface_la_LIBADD   = libtfbayes-dpm.la
dpm_gaussian_interface_la_LIBADD  += $(LIB_BOOST_PYTHON)
dpm_gaussian_interface_la_LDFLAGS  = -avoid-version -module

dpm_tfbs_interface_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
dpm_tfbs_interface_la_CXXFLAGS =
dpm_tfbs_interface_la_SOURCES  = \
	dpm-tfbs-interface.cc
dpm_tfbs_interface_la_LIBADD   = libtfbayes-dpm.la
dpm_tfbs_interface_la_LIBADD  += $(LIB_BOOST_PYTHON)
dpm_tfbs_interface_la_LDFLAGS  = -avoid-version -module

## extra dist
EXTRA_DIST = dpm-tfbs-demo-1.cfg dpm-tfbs-demo-1.fa.in dpm-tfbs-demo-2.cfg dpm-tfbs-demo-2.fa.in

## clean python files
clean-local:
	$(RM) *.pyc
