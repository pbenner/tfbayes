## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/include
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

bin_PROGRAMS = tfbayes-generate-alignment tfbayes-treespace-optimize tfbayes-treespace-mean tfbayes-treespace-histogram

tfbayes_generate_alignment_SOURCES = tfbayes-generate-alignment.cc
tfbayes_generate_alignment_LDADD   = ../phylotree/libtfbayes-phylotree.la

tfbayes_treespace_optimize_SOURCES = tfbayes-treespace-optimize.cc
tfbayes_treespace_optimize_LDADD   = -lpthread
tfbayes_treespace_optimize_LDADD  += ../phylotree/libtfbayes-phylotree.la
tfbayes_treespace_optimize_LDADD  += ../fasta/libtfbayes-fasta.la

tfbayes_treespace_mean_SOURCES = tfbayes-treespace-mean.cc
tfbayes_treespace_mean_LDADD   = ../phylotree/libtfbayes-phylotree.la

tfbayes_treespace_histogram_SOURCES = tfbayes-treespace-histogram.cc
tfbayes_treespace_histogram_LDADD   = ../phylotree/libtfbayes-phylotree.la

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/dpm

## scripts
dist_bin_SCRIPTS = \
	tfbayes-align \
	tfbayes-approximate \
	tfbayes-bed \
	tfbayes-expectation \
	tfbayes-cluster \
	tfbayes-generate \
	tfbayes-mean \
	tfbayes-motifformat \
	tfbayes-fasta-to-maf \
	tfbayes-maf-to-fasta \
	tfbayes-optimize \
	tfbayes-plot \
	tfbayes-partition \
	tfbayes-preprocess-alignment \
	tfbayes-sampler \
	tfbayes-scan

## clean python files
clean-local:
	$(RM) *.pyc
