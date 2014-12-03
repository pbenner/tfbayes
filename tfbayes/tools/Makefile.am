## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)
AM_CXXFLAGS = $(CXXFLAGS_NORTTI)

bin_PROGRAMS = tfbayes-generate-alignment tfbayes-treespace-sampler tfbayes-hmm tfbayes-treespace-estimate tfbayes-treespace-histogram tfbayes-fifo

tfbayes_fifo_SOURCES = tfbayes-fifo.cc

tfbayes_generate_alignment_SOURCES = tfbayes-generate-alignment.cc
tfbayes_generate_alignment_LDADD   = $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la

tfbayes_treespace_sampler_SOURCES  = tfbayes-treespace-sampler.cc
tfbayes_treespace_sampler_LDADD    = $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la
tfbayes_treespace_sampler_LDADD   += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
tfbayes_treespace_sampler_LDADD   += $(top_builddir)/tfbayes/fastarithmetics/libtfbayes-fastarithmetics.la

tfbayes_treespace_estimate_SOURCES = tfbayes-treespace-estimate.cc
tfbayes_treespace_estimate_LDADD   = $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la

tfbayes_treespace_histogram_SOURCES = tfbayes-treespace-histogram.cc
tfbayes_treespace_histogram_LDADD   = $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la

tfbayes_hmm_SOURCES = tfbayes-hmm.cc
tfbayes_hmm_LDADD   = -lm
tfbayes_hmm_LDADD  += $(top_builddir)/tfbayes/uipac/libtfbayes-uipac.la
tfbayes_hmm_LDADD  += $(top_builddir)/tfbayes/phylotree/libtfbayes-phylotree.la
tfbayes_hmm_LDADD  += $(top_builddir)/tfbayes/fasta/libtfbayes-fasta.la
tfbayes_hmm_LDADD  += $(BOOST_SYSTEM_LIB)

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/dpm

## scripts
dist_bin_SCRIPTS = \
	tfbayes-align \
	tfbayes-alignment-likelihood \
	tfbayes-alignment-shuffle \
	tfbayes-approximate \
	tfbayes-bed \
	tfbayes-estimate \
	tfbayes-expectation \
	tfbayes-export-results \
	tfbayes-motifformat \
	tfbayes-fasta-to-maf \
	tfbayes-locate \
	tfbayes-maf-to-fasta \
	tfbayes-plot \
	tfbayes-partition \
	tfbayes-preprocess-alignment \
	tfbayes-sampler \
	tfbayes-treespace-distance \
	tfbayes-treespace-geodesic \
	tfbayes-scan

## clean python files
clean-local:
	$(RM) *.pyc
