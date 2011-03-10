# Copyright (C) 2004, 2005, 2006, 2007 Philipp Benner <pbenner@uos.de>
#
# This makefile is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This makefile is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this makefile; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# There are quite some problems one has to deal with when using LaTeX and
# gnu make. First, there is no easy way to calculate dependencies for
# tex files. Second, several programs generate pdf output, hence there
# is no distinct target rule for pdf files. What one needs is a make
# system that knows how to distinguish between e.g. a pdf target generated
# from a dot graph or a tex file. Furthermore, the make system should
# know which files depend on each other according to their file names,
# to avoid having to parse tex files for dependencies. Another issue
# is that some files might exist or not. E.g. a tex document might have
# a bibliography based on bibtex and if that is the case the make
# system should recognize the corresponding bibtex file as dependency.

SHELL		= /bin/sh -e

LATEX		= latex -halt-on-error -shell-escape
TEXI2DVI	= texi2dvi -p
MAKEINDEX	= makeindex
BIBTEX		= bibtex

TEXINPUTS	= .:../share:../../share:

TEX_OUTPUT	= tfbs.pdf

## default build rule
all:
	@echo "Usage: make <pdf-output>"
	@echo "       make html"

## latex&co build rules
$(TEX_OUTPUT): $(wildcard $(basename $(TEX_OUTPUT))*.texinput)
$(TEX_OUTPUT): $(wildcard *.bib)
$(TEX_OUTPUT): %.pdf: %.tex
	TEXINPUTS=$(TEXINPUTS) $(LATEX) $(LATEXFLAGS) $*
#	$(BIBTEX) $*
#	TEXINPUTS=$(TEXINPUTS) $(LATEX) $(LATEXFLAGS) $*
	dvipdf $*.dvi

.SUFFIXES: .nlo .nls
.nlo.nls:
	$(MAKEINDEX) $< -s nomencl.ist -o $@
#	$(MAKEINDEX) $< -s nomentbl.ist -o $@

## clean rules
clean:
	$(RM) -r $(TEX_OUTPUT:.pdf=)
	$(RM) -r input/auto
	$(RM) -r auto
	$(RM) tmp.inputs
	$(RM) *.table
	$(RM) *.gnuplot
	$(RM) *.aux
	$(RM) *.blg
	$(RM) *.bbl
	$(RM) *.dvi
	$(RM) *.log
	$(RM) *.ps
	$(RM) *.toc
	$(RM) *.out
	$(RM) *.nav
	$(RM) *.snm
	$(RM) *.nlo
	$(RM) *.nls
	$(RM) *.ilg
	$(RM) *.idx
	$(RM) *.ind
	$(RM) *.bm

distclean: clean
	$(RM) $(TEX_OUTPUT)

.PHONY: all clean distclean html
