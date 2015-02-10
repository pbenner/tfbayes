## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

noinst_PROGRAMS = probability-test distribution-test

probability_test_SOURCES = probability-test.cc
distribution_test_SOURCES = distribution-test.cc

## headers
noinst_HEADERS = \
	abysmal-stack.hh \
	boost-random-shuffle.hh \
	clonable.hh \
	debug.hh \
	default-operator.hh \
	distribution.hh \
	histogram.hh \
	linalg.hh \
	logarithmetic.hh \
	named-ptr.hh \
	observable.hh \
	permutation.hh \
	polynomial.hh \
	probability.hh \
	progress.hh \
	random.hh \
	statistics.hh \
	strtools.hh \
	summation.hh \
	terminal-codes.hh \
	thread-pool.hh
