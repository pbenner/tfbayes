## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)

pkgpython_LTLIBRARIES = datatypes.la

datatypes_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
datatypes_la_SOURCES  = \
	datatypes.cc
datatypes_la_LIBADD   = $(LIB_BOOST_PYTHON)
datatypes_la_LDFLAGS  = -avoid-version -module

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/interface

## compile python files
pkgpython_PYTHON = __init__.py

## clean python files
clean-local:
	$(RM) *.pyc
