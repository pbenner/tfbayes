## Process this file with automake to produce Makefile.in

## set python package name
pkgpythondir = $(pyexecdir)/tfbayes/dpm

## scripts
dist_bin_SCRIPTS = \
	dpm-gaussian-demo \
	dpm-tfbs-align \
	dpm-tfbs-demo \
	dpm-tfbs-graph \
	dpm-tfbs-generate \
	dpm-tfbs-motifformat \
	dpm-tfbs-plot \
	dpm-tfbs-partition \
	dpm-tfbs-sampler \
	dpm-tfbs-scan \
	dpm-tfbs-tools \
	dpm-tfbs-weblogo \
	macs2multiz

## extra dist
EXTRA_DIST = dpm-tfbs-demo-1.cfg dpm-tfbs-demo-1.fa dpm-tfbs-demo-2.cfg dpm-tfbs-demo-2.fa

## clean python files
clean-local:
	$(RM) *.pyc
