## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = libexception include humandb

EXTRA_DIST = hdb-load-seq.in hdb-read-seq.in

dist_bin_SCRIPTS = hdb-load-seq hdb-read-seq

do_subst = sed -e 's,[@]pythondir[@],$(pythondir),g'

hdb-load-seq: hdb-load-seq.in Makefile
	$(do_subst) < $(srcdir)/hdb-load-seq.in > hdb-load-seq
	chmod +x hdb-load-seq

hdb-read-seq: hdb-read-seq.in Makefile
	$(do_subst) < $(srcdir)/hdb-read-seq.in > hdb-read-seq
	chmod +x hdb-read-seq
