## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

bin_PROGRAMS = polynomial-test likelihood-test marginal-likelihood-test

polynomial_test_SOURCES = polynomial-test.cc
polynomial_test_LDADD   = -lm

likelihood_test_SOURCES = likelihood-test.cc
likelihood_test_LDADD   = -lm
likelihood_test_LDADD  += libphylotree.la

marginal_likelihood_test_SOURCES = marginal-likelihood-test.cc
marginal_likelihood_test_LDADD   = -lm
marginal_likelihood_test_LDADD  += libphylotree.la

noinst_LTLIBRARIES = libphylotree.la

libphylotree_la_SOURCES = phylotree.cc phylotree.hh likelihood.cc likelihood.hh marginal-likelihood.cc marginal-likelihood.hh utility.cc utility.hh
libphylotree_la_LIBADD  = -lm
