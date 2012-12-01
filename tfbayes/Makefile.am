## Process this file with automake to produce Makefile.in

SUBDIRS  = fasta
#SUBDISS += humandb
SUBDIRS += dpm
SUBDIRS += motif
SUBDIRS += phylotree
SUBDIRS += uipac
SUBDIRS += tools

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes

## compile python files
pkgpython_PYTHON = __init__.py

## clean python files
clean-local:
	$(RM) *.pyc
