## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)

## install headers
otherincludedir = $(includedir)/tfbayes/interface
otherinclude_HEADERS = common.hh

## library
lib_LTLIBRARIES = libtfbayes-interface.la

libtfbayes_interface_la_SOURCES = \
	common.cc
libtfbayes_interface_la_LDFLAGS = -no-undefined -version-info 0:0:0

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/interface

## compile python files
pkgpython_PYTHON = __init__.py datatypes.py tools.py

## clean python files
clean-local:
	$(RM) *.pyc
