## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

noinst_PROGRAMS = entropy-test entropy-approximation

entropy_test_SOURCES = entropy-test.cc

entropy_approximation_SOURCES = entropy-approximation.cc
entropy_approximation_LDADD = $(BOOST_SERIALIZATION_LIB)

## headers
noinst_HEADERS = \
	entropy.hh
