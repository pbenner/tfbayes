## Process this file with automake to produce Makefile.in
## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

noinst_PROGRAMS = distribution-test

distribution_test_SOURCES = distribution-test.cc

## fg library
# lib_LTLIBRARIES = libtfbayes-fg.la
# libtfbayes_fg_la_SOURCES = \
# 	distribution.hh \
# 	distribution.cc \
# 	factor-graph.cc \
# 	factor-graph.hh
# libtfbayes_fg_la_LIBADD   = $(LIB_PTHREAD)
# libtfbayes_fg_la_LIBADD  += $(BOOST_REGEX_LIB)
# libtfbayes_fg_la_LIBADD  += $(BOOST_SYSTEM_LIB)
# libtfbayes_fg_la_LIBADD  += $(BOOST_THREAD_LIB)
# libtfbayes_fg_la_LDFLAGS  = -no-undefined
# libtfbayes_fg_la_LDFLAGS += -version-info 0:0:0

## python interface
pkgpython_LTLIBRARIES = interface.la

interface_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
interface_la_CXXFLAGS = 
interface_la_SOURCES  = \
	interface.cc
interface_la_LIBADD   = $(BOOST_PYTHON_LIB)
interface_la_LIBADD  += -lm
interface_la_LDFLAGS  = -avoid-version
interface_la_LDFLAGS += -module

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/fg

## compile python files
pkgpython_PYTHON = __init__.py

## clean python files
clean-local:
	$(RM) *.pyc
