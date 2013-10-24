## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)

## install headers
otherincludedir = $(includedir)/tfbayes/interface
otherinclude_HEADERS = exceptions.hh utility.hh

## interface
pkgpython_LTLIBRARIES = datatypes.la

datatypes_la_CPPFLAGS  = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
datatypes_la_SOURCES   = \
	datatypes.cc
datatypes_la_LIBADD    = $(BOOST_PYTHON_LIB)
datatypes_la_LDFLAGS   = -avoid-version
datatypes_la_LDFLAGS  += -module

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/interface

## compile python files
pkgpython_PYTHON = __init__.py sigint.py

## clean python files
clean-local:
	$(RM) *.pyc
