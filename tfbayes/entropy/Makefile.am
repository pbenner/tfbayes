## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

#noinst_PROGRAMS = entropy-test entropy-approximation entropy-approximation-test entropy-distribution-test entropy-multinomial-distribution-test

entropy_test_SOURCES = entropy-test.cc

entropy_approximation_SOURCES = \
	entropy-approximation-mixture.hh \
	entropy-approximation-recursive.hh \
	entropy-approximation.cc
entropy_approximation_LDADD = $(BOOST_SERIALIZATION_LIB)

entropy_approximation_test_SOURCES = entropy-approximation-test.cc
entropy_approximation_test_LDADD = $(BOOST_SERIALIZATION_LIB)

entropy_distribution_test_SOURCES = entropy-distribution-test.cc
entropy_distribution_test_LDADD = $(BOOST_SERIALIZATION_LIB)

entropy_multinomial_distribution_test_SOURCES = entropy-multinomial-distribution-test.cc
entropy_multinomial_distribution_test_LDADD = $(BOOST_SERIALIZATION_LIB)


## headers
noinst_HEADERS = \
	entropy-approximation.hh \
	entropy-approximation-2.hh \
	entropy-approximation-3.hh \
	entropy-approximation-4.hh \
	entropy-approximation-5.hh \
	entropy-approximation-6.hh \
	entropy-approximation-7.hh \
	entropy-approximation-8.hh \
	entropy-approximation-9.hh \
	entropy-approximation-10.hh \
	entropy-distribution.hh \
	entropy-multinomial-distribution.hh \
	entropy.hh \
	multinomial-distribution.hh
