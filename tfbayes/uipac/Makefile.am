## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

noinst_PROGRAMS = alphabet-test

alphabet_test_SOURCES = alphabet-test.cc
alphabet_test_LDADD   = libtfbayes-uipac.la

## library
lib_LTLIBRARIES = libtfbayes-uipac.la

libtfbayes_uipac_la_SOURCES = \
	alphabet.cc
	alphabet.hh
libtfbayes_uipac_la_LDFLAGS = -no-undefined -version-info 0:0:0

## python interface
pkgpython_LTLIBRARIES = interface.la

interface_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
interface_la_CXXFLAGS = 
interface_la_SOURCES  = \
	interface.cc
interface_la_LIBADD   = $(LIB_BOOST_PYTHON)
interface_la_LIBADD  += libtfbayes-uipac.la
interface_la_LDFLAGS  = -avoid-version -module

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/uipac

## compile python files
pkgpython_PYTHON = __init__.py alphabet.py

## clean python files
clean-local:
	$(RM) *.pyc
