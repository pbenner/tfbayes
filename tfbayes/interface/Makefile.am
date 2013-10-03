## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)

## install headers
otherincludedir = $(includedir)/tfbayes/interface
otherinclude_HEADERS = datatypes.hh

## library
lib_LTLIBRARIES = libtfbayes-interface.la

pkgpython_LTLIBRARIES = linalg.la

linalg_la_CPPFLAGS = $(PYTHON_CPPFLAGS)
linalg_la_SOURCES  = \
	linalg.cc
linalg_la_LIBADD   = $(LIB_BOOST_PYTHON)
linalg_la_LDFLAGS  = -avoid-version -module

libtfbayes_interface_la_SOURCES = \
	datatypes.cc
libtfbayes_interface_la_LDFLAGS = -no-undefined -version-info 0:0:0

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/interface

## compile python files
pkgpython_PYTHON = __init__.py datatypes.py loadlibrary.py

## clean python files
clean-local:
	$(RM) *.pyc
