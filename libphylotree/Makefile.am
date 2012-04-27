## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## flex/bison
LEX  = flex
YACC = bison -y
BUILT_SOURCES = phylotree-parser.h
AM_YFLAGS = -d
AM_LFLAGS = -o$(LEX_OUTPUT_ROOT).c

bin_PROGRAMS = polynomial-test phylotree-expand-test phylotree-parser-test likelihood-test likelihood-prime-test marginal-likelihood-test

polynomial_test_SOURCES = polynomial-test.cc
polynomial_test_LDADD   = -lm

phylotree_expand_test_SOURCES = phylotree-expand-test.cc
phylotree_expand_test_LDADD   = -lm
phylotree_expand_test_LDADD  += libphylotree.la

phylotree_parser_test_SOURCES = phylotree-parser-test.cc
phylotree_parser_test_LDADD   = -lm
phylotree_parser_test_LDADD  += libphylotree.la

likelihood_test_SOURCES = likelihood-test.cc
likelihood_test_LDADD   = -lm
likelihood_test_LDADD  += libphylotree.la

likelihood_prime_test_SOURCES = likelihood-prime-test.cc
likelihood_prime_test_LDADD   = -lm
likelihood_prime_test_LDADD  += libphylotree.la

marginal_likelihood_test_SOURCES = marginal-likelihood-test.cc
marginal_likelihood_test_LDADD   = -lm
marginal_likelihood_test_LDADD  += libphylotree.la

noinst_LTLIBRARIES = libphylotree.la

libphylotree_la_SOURCES = \
	nodeset.cc \
	nodeset.hh \
	phylotree.cc \
	phylotree.hh \
	phylotree-expand.hh \
	phylotree-polynomial.hh \
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
	utility.cc \
	utility.hh
libphylotree_la_LIBADD  = -lm
