## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/include
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

## library
noinst_LTLIBRARIES = libtfbayes-exception.la

libtfbayes_exception_la_SOURCES = exception.c exception.h

## install headers
otherincludedir = $(includedir)/tfbayes/exception
otherinclude_HEADERS = \
	exception.h
