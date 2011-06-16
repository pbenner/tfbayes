## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = libexception include humandb dpm

EXTRA_DIST = hdb-load-seq.in hdb-load-chrom.in hdb-read-seq.in hdb-read-chrom.in hdb-search-seq.in hdb-statistics-chrom.in hdb-statistics-seq.in

dist_bin_SCRIPTS = hdb-load-seq hdb-load-chrom hdb-read-seq hdb-read-chrom hdb-search-seq hdb-statistics-chrom hdb-statistics-seq

do_subst = sed -e 's,[@]pythondir[@],$(pythondir),g'

hdb-load-seq: hdb-load-seq.in Makefile
	$(do_subst) < $(srcdir)/hdb-load-seq.in > hdb-load-seq
	chmod +x hdb-load-seq

hdb-load-chrom: hdb-load-chrom.in Makefile
	$(do_subst) < $(srcdir)/hdb-load-chrom.in > hdb-load-chrom
	chmod +x hdb-load-chrom

hdb-read-seq: hdb-read-seq.in Makefile
	$(do_subst) < $(srcdir)/hdb-read-seq.in > hdb-read-seq
	chmod +x hdb-read-seq

hdb-read-chrom: hdb-read-chrom.in Makefile
	$(do_subst) < $(srcdir)/hdb-read-chrom.in > hdb-read-chrom
	chmod +x hdb-read-chrom

hdb-search-seq: hdb-search-seq.in Makefile
	$(do_subst) < $(srcdir)/hdb-search-seq.in > hdb-search-seq
	chmod +x hdb-search-seq

hdb-statistics-chrom: hdb-statistics-chrom.in Makefile
	$(do_subst) < $(srcdir)/hdb-statistics-chrom.in > hdb-statistics-chrom
	chmod +x hdb-statistics-chrom

hdb-statistics-seq: hdb-statistics-seq.in Makefile
	$(do_subst) < $(srcdir)/hdb-statistics-seq.in > hdb-statistics-seq
	chmod +x hdb-statistics-seq
