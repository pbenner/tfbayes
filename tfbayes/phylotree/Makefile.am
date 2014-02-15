## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

## flex/bison
LEX  = flex
YACC = bison -y
BUILT_SOURCES = parser.hh
AM_YFLAGS = -d
AM_LFLAGS = -o$(LEX_OUTPUT_ROOT).c

noinst_PROGRAMS = polynomial-test tree-reduction-test phylotree-test expected-terms parser-test approximation-test likelihood-test marginal-likelihood-test treespace-test treespace-distance

polynomial_test_SOURCES = polynomial-test.cc
polynomial_test_LDADD   = -lm

tree_reduction_test_SOURCES = tree-reduction-test.cc
tree_reduction_test_LDADD   = -lm
tree_reduction_test_LDADD  += libtfbayes-phylotree.la

phylotree_test_SOURCES = phylotree-test.cc
phylotree_test_LDADD   = -lm
phylotree_test_LDADD  += libtfbayes-phylotree.la

expected_terms_SOURCES = expected-terms.cc
expected_terms_LDADD   = -lm
expected_terms_LDADD  += libtfbayes-phylotree.la

parser_test_SOURCES = parser-test.cc
parser_test_LDADD   = -lm
parser_test_LDADD  += libtfbayes-phylotree.la

approximation_test_SOURCES = approximation-test.cc
approximation_test_LDADD   = -lm
approximation_test_LDADD  += libtfbayes-phylotree.la

likelihood_test_SOURCES = likelihood-test.cc
likelihood_test_LDADD   = -lm
likelihood_test_LDADD  += libtfbayes-phylotree.la

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
	approximation.hh \
	generate-observations.hh \
	gradient.hh \
	gradient-ascent.hh \
	lexer.ll \
	model.hh \
	parser.hh \
	parser.yy \
	parsetree.cc \
	parsetree.hh \
	phylotree.cc \
	phylotree.hh \
	polynomial.hh \
	posterior.hh \
	marginal-likelihood.hh \
	utility.hh \
	sampler.hh \
	simple-polynomial.hh \
	treespace.cc \
	treespace.hh
libtfbayes_phylotree_la_LIBADD   = -lm
libtfbayes_phylotree_la_LIBADD  += $(LIB_GLPK)
libtfbayes_phylotree_la_LIBADD  += $(LIB_PTHREAD)
libtfbayes_phylotree_la_LIBADD  += $(BOOST_SYSTEM_LIB)
libtfbayes_phylotree_la_LIBADD  += $(BOOST_THREAD_LIB)
libtfbayes_phylotree_la_LIBADD  += $(top_builddir)/tfbayes/uipac/libtfbayes-uipac.la
libtfbayes_phylotree_la_LDFLAGS  = -no-undefined 
libtfbayes_phylotree_la_LDFLAGS += -version-info 0:0:0

## other files that should be distributed
EXTRA_DIST = phylotree.nb

pkgpython_LTLIBRARIES = interface.la

interface_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
interface_la_CXXFLAGS = 
interface_la_SOURCES  = \
	interface.cc
interface_la_LIBADD   = $(BOOST_PYTHON_LIB)
interface_la_LIBADD  += libtfbayes-phylotree.la
interface_la_LDFLAGS  = -avoid-version
interface_la_LDFLAGS += -module

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/phylotree

## compile python files
pkgpython_PYTHON = __init__.py

## clean python files
clean-local:
	$(RM) *.pyc
