## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = libexception include humandb dpm

EXTRA_DIST = hdb-load-upstream.in hdb-load-chrom.in hdb-read-upstream.in hdb-read-chrom.in hdb-search-upstream.in hdb-statistics-chrom.in hdb-statistics-upstream.in

dist_bin_SCRIPTS = hdb-load-upstream hdb-load-chrom hdb-read-upstream hdb-read-chrom hdb-search-upstream hdb-statistics-chrom hdb-statistics-upstream

do_subst = sed -e 's,[@]pythondir[@],$(pythondir),g'

hdb-load-upstream: hdb-load-upstream.in Makefile
	$(do_subst) < $(srcdir)/hdb-load-upstream.in > hdb-load-upstream
	chmod +x hdb-load-upstream

hdb-load-chrom: hdb-load-chrom.in Makefile
	$(do_subst) < $(srcdir)/hdb-load-chrom.in > hdb-load-chrom
	chmod +x hdb-load-chrom

hdb-read-upstream: hdb-read-upstream.in Makefile
	$(do_subst) < $(srcdir)/hdb-read-upstream.in > hdb-read-upstream
	chmod +x hdb-read-upstream

hdb-read-chrom: hdb-read-chrom.in Makefile
	$(do_subst) < $(srcdir)/hdb-read-chrom.in > hdb-read-chrom
	chmod +x hdb-read-chrom

hdb-search-upstream: hdb-search-upstream.in Makefile
	$(do_subst) < $(srcdir)/hdb-search-upstream.in > hdb-search-upstream
	chmod +x hdb-search-upstream

hdb-statistics-chrom: hdb-statistics-chrom.in Makefile
	$(do_subst) < $(srcdir)/hdb-statistics-chrom.in > hdb-statistics-chrom
	chmod +x hdb-statistics-chrom

hdb-statistics-upstream: hdb-statistics-upstream.in Makefile
	$(do_subst) < $(srcdir)/hdb-statistics-upstream.in > hdb-statistics-upstream
	chmod +x hdb-statistics-upstream
