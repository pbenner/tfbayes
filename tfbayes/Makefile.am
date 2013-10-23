## Process this file with automake to produce Makefile.in

SUBDIRS  = alignio
SUBDIRS += weblogo
SUBDIRS += corebio
SUBDIRS += interface
SUBDIRS += fastarithmetics
SUBDIRS += fasta
SUBDIRS += uipac
SUBDIRS += dpm
SUBDIRS += config
SUBDIRS += cluster
SUBDIRS += phylotree
SUBDIRS += alignment
SUBDIRS += utility
SUBDIRS += tools

otherincludedir = $(includedir)/tfbayes
otherinclude_HEADERS = \
	config.h

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes

## compile python files
pkgpython_PYTHON = __init__.py

## clean python files
clean-local:
	$(RM) *.pyc
