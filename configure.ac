[1]dnl -- configure.ac --
dnl
dnl Copyright (C) 2010 Philipp Benner
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

AC_INIT([TFBayes], [1.0], [philipp.benner@mis.mpg.de])
AC_CONFIG_AUX_DIR(auxtools)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([auxtools/compile])
AM_INIT_AUTOMAKE([gnu dist-bzip2])
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL
LT_INIT
AC_PROG_CC
AC_PROG_CXX
AC_C_CONST
AC_ISC_POSIX
AC_PROG_INSTALL
AM_PROG_CC_C_O
AC_PROG_LN_S

dnl ,---------------------------- 
dnl | GENERAL PACKAGE INFO
dnl `----------------------------

AC_SUBST([PACKAGE_SUMMARY], ["template project"])
AC_SUBST([PACKAGE_URL], ["http://www.homepage.ch"])
AC_SUBST([PACKAGE_LICENSE], ["GPL"])
AC_SUBST([PACKAGE_DESCRIPTION], ["This is a template project using autoconf, automake and autoheader."])

dnl ,---------------------------- 
dnl | MACRO DEFINITIONS
dnl `----------------------------

AC_DEFUN([AC_MY_ARG_WITH],
	[# create help entry and conditionals
		AC_ARG_WITH([$1],
			AC_HELP_STRING([--with-$1],
				[# help string
					build $1 (default is $2)]),
			ac_cv_use_$1=$withval,
			ac_cv_use_$1=$2)
		AC_CACHE_CHECK(whether to use $1, ac_cv_use_$1, ac_cv_use_$1=$2)])

dnl AC_NEW_SUBDIR(subdir_name, default_config)
AC_DEFUN([AC_NEW_SUBDIR],
        [AC_MY_ARG_WITH([$1], [$3])
	AM_CONDITIONAL([COND_$2], [test "$ac_cv_use_$1" = yes])])

dnl ,---------------------------- 
dnl | LIBRARIES
dnl `----------------------------

dnl HP-UX provides a library called libhplx which provides linux apis.
AC_ARG_VAR(HPLX_LIBRARY_PATH, [Path to the hplx library directory (only for hp-ux systems).])
AC_ARG_VAR(HPLX_INCLUDE_PATH, [Path to the hplx include directory (only for hp-ux systems).])
AC_CHECK_LIB([hplx], [hstrerror],
	[# Add additional include an library directory.
	LDFLAGS="$LDFLAGS -L${HPLX_LIBRARY_PATH:-/usr/local/hplx/lib}"
	CFLAGS="$CFLAGS -I${HPLX_INCLUDE_PATH:-/usr/local/hplx/include}"],,
	[# Additional libraries:
	-lsec])

AC_CHECK_LIB([ucb], [sprintf])
dnl Tru64 distributes vasprintf() in -liberty
AC_CHECK_LIB([iberty], [vasprintf])
AC_CHECK_LIB([m],[cos],,
  AC_MSG_ERROR([this library is required]))
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,
  AC_MSG_ERROR([this library is required]))
AC_CHECK_LIB([gsl], [gsl_blas_dgemm],,
  AC_MSG_ERROR([this library is required]))
AC_CHECK_LIB([db], [db_env_set_func_exists], LIB_DB="-ldb",
  AC_MSG_ERROR([this library is required]))
AC_CHECK_LIB([pthread], [pthread_create], LIB_PTHREAD="-lpthread",
  AC_MSG_ERROR([this library is required]))

AC_SUBST(LIB_DB)
AC_SUBST(LIB_PTHREAD)

dnl ,---------------------------- 
dnl | GETOPT_LONG
dnl `----------------------------

AC_DEFUN([FUNC_GETOPT_LONG],
	[AC_PREREQ(2.49)dnl
		# clean out junk possibly left behind by a previous configuration
		rm -f include/getopt.h
		# Check for getopt_long support
		AC_CHECK_HEADERS([getopt.h])
		AC_CHECK_FUNCS([getopt_long],,
			[# FreeBSD has a gnugetopt library for this
				AC_CHECK_LIB([gnugetopt],[getopt_long],[AC_DEFINE([HAVE_GETOPT_LONG])],
					[# use the GNU replacement
						AC_LIBOBJ(getopt)
						AC_LIBOBJ(getopt1)
						AC_CONFIG_LINKS([include/getopt.h:include/gnugetopt.h])])])])

FUNC_GETOPT_LONG

dnl ,---------------------------- 
dnl | DEFINES
dnl `----------------------------

AH_VERBATIM([_GNU_SOURCE],
[/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif])

AC_DEFINE_UNQUOTED(TARGET, "$target",
	[System type this program will compile code for (only needed for compilers).])
AC_DEFINE_UNQUOTED(TARGET_CPU, "$target_cpu",
	[System cpu this program will compile code for (only needed for compilers).])
AC_DEFINE_UNQUOTED(TARGET_OS, "$target_os",
	[Operatins system this program will compile code for (only needed for compilers).])
AC_DEFINE_UNQUOTED(TARGET_VENDOR, "$target_vendor",
	[System vendor type this program will compile code for (only needed for compilers).])

AC_DEFINE_UNQUOTED(BUILD, "$build",
	[System type this program was built on.])
AC_DEFINE_UNQUOTED(BUILD_CPU, "$build_cpu",
	[System cpu this program was built on.])
AC_DEFINE_UNQUOTED(BUILD_OS, "$build_os",
	[Operatins system this program was built on.])
AC_DEFINE_UNQUOTED(BUILD_VENDOR, "$build_vendor",
	[System vendor type this program was built on.])

AC_DEFINE_UNQUOTED(HOST, "$host",
	[System type this program will run on.])
AC_DEFINE_UNQUOTED(HOST_CPU, "$host_cpu",
	[System cpu this program will run on.])
AC_DEFINE_UNQUOTED(HOST_OS, "$host_os",
	[Operatins system this program will run on.])
AC_DEFINE_UNQUOTED(HOST_VENDOR, "$host_vendor",
	[System vendor type this program will run on.])

dnl ,---------------------------- 
dnl | HEADERS
dnl `----------------------------

HEADERS_TO_CHECK="stdarg.h dlfcn.h dirent.h netinet/in.h sys/ioctl.h
									arpa/inet.h netdb.h syslog.h signal.h pthread.h"

AC_CHECK_HEADERS($HEADERS_TO_CHECK,,
	AC_MSG_ERROR([this header is required]))

AC_CHECK_HEADER(net/if.h,,
	AC_MSG_ERROR([this header is required]),
	[[#include <sys/types.h>
		#include <sys/socket.h>]] )

AC_CHECK_HEADER([error.h],
	AC_DEFINE(HAVE_ERROR_H, [], [Needed on OpenBSD systems.]))

dnl ,---------------------------- 
dnl | FUNCTIONS
dnl `----------------------------

AC_CHECK_FUNCS(dprintf)
AC_CHECK_FUNCS(vsyslog)

dnl ,---------------------------- 
dnl | PRECISION
dnl `----------------------------

AC_ARG_ENABLE([longdouble],
	[AS_HELP_STRING([--enable-longdouble],
								[use long double for probabilities @<:@default=yes@:>@])],
	[],
	[enable_longdouble=yes])
AS_IF([test "x$enable_longdouble" = "xyes"],
	[AC_CHECK_FUNCS([expl logl log1pl], [], [enable_longdouble=no])], [])
AS_IF([test "x$enable_longdouble" = "xyes"],
	[AC_DEFINE([PROB_T], [long double], [Use long double to store log probabilities.])],
	[AC_DEFINE([PROB_T], [double], [Use double to store log probabilities.])
	 AC_DEFINE([expl], [exp], [Replace expl with exp])
	 AC_DEFINE([logl], [log], [Replace logl with log])
	 AC_DEFINE([log1pl], [log1p], [Replace log1pl with log1p])
	 ])

dnl ,---------------------------- 
dnl | GCC_OPTIONS
dnl `----------------------------

CFLAGS=""
CXXFLAGS=""
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
								[turn on debugging @<:@default=no@:>@])],
	[],
	[enable_debug=no])
AS_IF([test "x$enable_debug" = "xyes"],
  [AX_CFLAGS_GCC_OPTION([-O0])
   AX_CFLAGS_GCC_OPTION([-g])
   AX_CFLAGS_GCC_OPTION([-ggdb])
   AX_CFLAGS_GCC_OPTION([-pg])
	 dnl
	 AX_CXXFLAGS_GCC_OPTION([-O0])
	 AX_CXXFLAGS_GCC_OPTION([-g])
   AX_CXXFLAGS_GCC_OPTION([-ggdb])
   AX_CXXFLAGS_GCC_OPTION([-pg])
	 CFLAGS="$CFLAGS -DDEBUG"
	 CXXFLAGS="$CXXFLAGS -DDEBUG"],
  [AX_CFLAGS_GCC_OPTION([-O4])
   AX_CFLAGS_GCC_OPTION([-ffast-math])
	 dnl
	 AX_CXXFLAGS_GCC_OPTION([-O4])
   AX_CXXFLAGS_GCC_OPTION([-ffast-math])])

AX_CFLAGS_GCC_OPTION([-Wall])
AX_CFLAGS_GCC_OPTION([-Wwrite-strings])
AX_CFLAGS_GCC_OPTION([-Wstrict-prototypes])
AX_CFLAGS_GCC_OPTION([-Wno-trigraphs])
AX_CFLAGS_GCC_OPTION([-Wtrampolines])
AX_CFLAGS_GCC_OPTION([-fPIC])
AX_CFLAGS_GCC_OPTION([-fno-nested-functions])

AX_CXXFLAGS_GCC_OPTION([-Wall])
AX_CXXFLAGS_GCC_OPTION([-Wwrite-strings])
AX_CXXFLAGS_GCC_OPTION([-Wstrict-prototypes])
AX_CXXFLAGS_GCC_OPTION([-Wno-trigraphs])
AX_CXXFLAGS_GCC_OPTION([-Wtrampolines])
AX_CXXFLAGS_GCC_OPTION([-fPIC])
AX_CXXFLAGS_GCC_OPTION([-fno-nested-functions])

dnl ,---------------------------- 
dnl | SUBDIR CONFIGURATION
dnl `----------------------------

#AC_NEW_SUBDIR([doc], [DOC], ["no"])

dnl ,---------------------------- 
dnl | PERL
dnl `----------------------------

dnl AC_PATH_PROG(PERL, perl, "no")
dnl AM_CONDITIONAL([COND_PERL], [test "$PERL" != no])

dnl ,---------------------------- 
dnl | PYTHON
dnl `----------------------------

AM_PATH_PYTHON(2.5)

dnl ,---------------------------- 
dnl | CYGWIN
dnl `----------------------------

dnl disable library dependency check on cygwin
dnl because of problems with -liberty
AS_IF([test "x$host_os" = "xcygwin"],
  [deplibs_check_method="pass_all"], [])

dnl ,---------------------------- 
dnl | AC_OUTPUT
dnl `----------------------------

AC_CONFIG_HEADER([config.h:config.h.in])
AC_CONFIG_FILES([
	Makefile
	version.h
	libexception/Makefile
	include/Makefile
	humandb/Makefile
	dpm/Makefile
 	])

AC_OUTPUT

dnl ,---------------------------- 
dnl | Print a summary
dnl `----------------------------

AC_MSG_NOTICE_BOX([Configure summary

System type:              ${host}
Build with debug symbols: ${enable_debug}
Use long double:          ${enable_longdouble}
])

dnl Local Variables:
dnl tab-width: 2
dnl End:
