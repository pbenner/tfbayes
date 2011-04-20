## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## Berkeley DB for the human genome data
lib_LTLIBRARIES  = libhumandb.la
libhumandb_la_SOURCES = humandb.c humandb.h interface.c interface.h
libhumandb_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libhumandb_la_LDFLAGS = -no-undefined -version-info 0:0:0

## clean python files
clean-local:
	$(RM) *.pyc
