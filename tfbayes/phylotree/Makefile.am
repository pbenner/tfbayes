## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)

## flex/bison
LEX  = flex
YACC = bison -y
BUILT_SOURCES = phylotree-parser.h
AM_YFLAGS = -d
AM_LFLAGS = -o$(LEX_OUTPUT_ROOT).c

bin_PROGRAMS = polynomial-test phylotree-test phylotree-expand-test phylotree-parser-test phylotree-approximation-test likelihood-test likelihood-prime-test marginal-likelihood-test treespace-test treespace-distance

polynomial_test_SOURCES = polynomial-test.cc
polynomial_test_LDADD   = -lm

phylotree_test_SOURCES = phylotree-test.cc
phylotree_test_LDADD   = -lm
phylotree_test_LDADD  += libtfbayes-phylotree.la

phylotree_expand_test_SOURCES = phylotree-expand-test.cc
phylotree_expand_test_LDADD   = -lm
phylotree_expand_test_LDADD  += libtfbayes-phylotree.la

# phylotree_expected_terms_SOURCES = phylotree-expected-terms.cc
# phylotree_expected_terms_LDADD   = -lm
# phylotree_expected_terms_LDADD  += libtfbayes-phylotree.la

phylotree_parser_test_SOURCES = phylotree-parser-test.cc
phylotree_parser_test_LDADD   = -lm
phylotree_parser_test_LDADD  += libtfbayes-phylotree.la

phylotree_approximation_test_SOURCES = phylotree-approximation-test.cc
phylotree_approximation_test_LDADD   = -lm
phylotree_approximation_test_LDADD  += libtfbayes-phylotree.la

likelihood_test_SOURCES = likelihood-test.cc
likelihood_test_LDADD   = -lm
likelihood_test_LDADD  += libtfbayes-phylotree.la

likelihood_prime_test_SOURCES = likelihood-prime-test.cc
likelihood_prime_test_LDADD   = -lm
likelihood_prime_test_LDADD  += libtfbayes-phylotree.la

marginal_likelihood_test_SOURCES = marginal-likelihood-test.cc
marginal_likelihood_test_LDADD   = -lm
marginal_likelihood_test_LDADD  += libtfbayes-phylotree.la

treespace_test_SOURCES = treespace-test.cc
treespace_test_LDADD   = -lm
treespace_test_LDADD  += libtfbayes-phylotree.la

treespace_distance_SOURCES = treespace-distance.cc
treespace_distance_LDADD   = -lm
treespace_distance_LDADD  += libtfbayes-phylotree.la

lib_LTLIBRARIES = libtfbayes-phylotree.la

libtfbayes_phylotree_la_SOURCES = \
	interface.hh \
	interface.cc \
	leafset.cc \
	leafset.hh \
	phylotree.cc \
	phylotree.hh \
	phylotree-approximation.hh \
	phylotree-expand.hh \
	phylotree-gradient.hh \
	phylotree-gradient-ascent.hh \
	phylotree-polynomial.hh \
	phylotree-sampler.hh \
	phylotree-lexer.ll \
	phylotree-parser.hh \
	phylotree-parser.yy \
	phylotree-parsetree.cc \
	phylotree-parsetree.hh \
	phylotree-simplify.cc \
	phylotree-simplify.hh \
	posterior.hh \
	incomplete-expression.cc \
	incomplete-expression.hh \
	marginal-likelihood.hh \
	statistics.hh \
	utility.cc \
	utility.hh \
	treespace.cc \
	treespace.hh
libtfbayes_phylotree_la_LIBADD  = -lm
libtfbayes_phylotree_la_LIBADD += $(LIB_GLPK)
libtfbayes_phylotree_la_LIBADD += $(top_builddir)/tfbayes/exception/libtfbayes-exception.la
libtfbayes_phylotree_la_LIBADD += $(top_builddir)/tfbayes/interface/libtfbayes-interface.la
libtfbayes_phylotree_la_LDFLAGS = -no-undefined -version-info 0:0:0

otherincludedir = $(includedir)/tfbayes/phylotree
otherinclude_HEADERS = \
	interface.hh \
	leafset.hh \
	phylotree.hh \
	phylotree-approximation.hh \
	phylotree-expand.hh \
	phylotree-gradient.hh \
	phylotree-gradient-ascent.hh \
	phylotree-polynomial.hh \
	phylotree-sampler.hh \
	phylotree-parser.h \
	phylotree-parser.hh \
	phylotree-parsetree.hh \
	phylotree-simplify.hh \
	posterior.hh \
	incomplete-expression.hh \
	marginal-likelihood.hh \
	statistics.hh \
	utility.hh \
	treespace.hh


## other files that should be distributed
EXTRA_DIST = phylotree-polynomial.nb

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/phylotree

## compile python files
pkgpython_PYTHON = __init__.py interface.py

## clean python files
clean-local:
	$(RM) *.pyc
