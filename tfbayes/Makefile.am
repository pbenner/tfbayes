## Process this file with automake to produce Makefile.in

SUBDIRS  = exception
SUBDIRS += interface
SUBDIRS += fastarithmetics
SUBDIRS += fasta
SUBDIRS += dpm
SUBDIRS += config
SUBDIRS += cluster
SUBDIRS += phylotree
SUBDIRS += alignment
SUBDIRS += uipac
SUBDIRS += tools
SUBDIRS += utility

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes

## compile python files
pkgpython_PYTHON = __init__.py

## clean python files
clean-local:
	$(RM) *.pyc
