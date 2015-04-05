## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

## flex/bison
LEX  = flex
YACC = bison -y
BUILT_SOURCES = partition-parser.hh
AM_YFLAGS = -d
AM_LFLAGS = -o$(LEX_OUTPUT_ROOT).c

# test programs
noinst_PROGRAMS = partition-parser-test

partition_parser_test_SOURCES = partition-parser-test.cc
partition_parser_test_LDADD   = -lm
partition_parser_test_LDADD  += libtfbayes-config.la
partition_parser_test_LDADD  += ../dpm/libtfbayes-dpm.la

# libraries
lib_LTLIBRARIES = libtfbayes-config.la
libtfbayes_config_la_SOURCES = \
	partition-code.hh \
	partition-code-requires.hh \
	partition-lexer.ll \
	partition-lexer.hh \
	partition-lexer-code.hh \
	partition-parser.yy

## interface
pkgpython_LTLIBRARIES = interface.la

interface_la_CPPFLAGS  = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
interface_la_SOURCES   = \
	interface.cc
interface_la_LIBADD    = libtfbayes-config.la
interface_la_LIBADD   += $(BOOST_PYTHON_LIB)
interface_la_LDFLAGS   = -avoid-version
interface_la_LDFLAGS  += -module

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/config

## compile python files
pkgpython_PYTHON = __init__.py results.py sampler.py parser.py tools.py

## clean python files
clean-local:
	$(RM) *.pyc
