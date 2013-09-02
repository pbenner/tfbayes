## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

noinst_PROGRAMS = dpm-tfbs-debug dpm-gaussian-debug

# mcm_test_SOURCES = mcm-test.cc
# mcm_test_LDADD   = libtfbayes-dpm.la
# mcm_test_LDADD  += -lm $(LIB_PTHREAD) $(LIB_BOOST_REGEX) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

dpm_gaussian_debug_SOURCES = dpm-gaussian-main.cc
dpm_gaussian_debug_LDADD   = $(top_builddir)/tfbayes/exception/libtfbayes-exception.la
dpm_gaussian_debug_LDADD  += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
dpm_gaussian_debug_LDADD  += libtfbayes-dpm.la
dpm_gaussian_debug_LDADD  += $(LIB_PTHREAD) $(LIB_BOOST_REGEX) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

dpm_tfbs_debug_SOURCES = dpm-tfbs-main.cc
dpm_tfbs_debug_LDADD   = $(top_builddir)/tfbayes/exception/libtfbayes-exception.la
dpm_tfbs_debug_LDADD  += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
dpm_tfbs_debug_LDADD  += libtfbayes-dpm.la
dpm_tfbs_debug_LDADD  += $(LIB_PTHREAD) $(LIB_BOOST_REGEX) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/dpm

## compile python files
pkgpython_PYTHON = __init__.py dpm_gaussian_interface.py dpm_tfbs_interface.py

## scripts
dist_bin_SCRIPTS = \
	dpm-gaussian-demo \
	dpm-tfbs-demo

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
	dpm-gaussian-interface.cc	\
	dpm-gaussian.cc			\
	dpm-gaussian.hh			\
	dpm-graph.hh			\
	dpm-tfbs.cc			\
	dpm-tfbs.hh			\
	dpm-tfbs-command.cc		\
	dpm-tfbs-command.hh		\
	dpm-tfbs-interface.cc		\
	dpm-tfbs-io.cc			\
	dpm-tfbs-io.hh			\
	dpm-tfbs-median.cc              \
	dpm-tfbs-median.hh              \
	dpm-tfbs-repl.cc		\
	dpm-tfbs-repl.hh		\
	dpm-tfbs-options.hh		\
	dpm-tfbs-prior.cc		\
	dpm-tfbs-prior.hh		\
	dpm-tfbs-sampler.cc		\
	dpm-tfbs-sampler.hh		\
	dpm-tfbs-state.cc		\
	dpm-tfbs-state.hh		\
	dpm-tfbs-test.cc		\
	dpm-partition.hh		\
	index.cc			\
	index.hh			\
	indexer.hh			\
	init.cc				\
	init.hh				\
	mixture-model.hh		\
	mixture-state.cc		\
	mixture-state.hh		\
	mixture-weights.hh		\
	observer.hh			\
	pmcmc.cc			\
	pmcmc.hh			\
	sampler.cc			\
	sampler.hh			\
	state.cc			\
	state.hh			\
	statistics.hh			\
	utility.hh
libtfbayes_dpm_la_LIBADD  = $(top_builddir)/tfbayes/exception/libtfbayes-exception.la
libtfbayes_dpm_la_LIBADD += $(top_builddir)/tfbayes/interface/libtfbayes-interface.la
libtfbayes_dpm_la_LIBADD += $(top_builddir)/tfbayes/fastarithmetics/libtfbayes-fastarithmetics.la
libtfbayes_dpm_la_LIBADD += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
libtfbayes_dpm_la_LIBADD += $(LIB_PTHREAD) $(LIB_BOOST_REGEX) $(LIB_BOOST_SYSTEM) $(LIB_BOOST_THREAD)
libtfbayes_dpm_la_LDFLAGS = -no-undefined -version-info 0:0:0

## extra dist
EXTRA_DIST = dpm-tfbs-demo-1.cfg dpm-tfbs-demo-1.fa.in dpm-tfbs-demo-2.cfg dpm-tfbs-demo-2.fa.in

## clean python files
clean-local:
	$(RM) *.pyc
