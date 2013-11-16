[1]dnl -- configure.ac --
dnl
dnl Copyright (C) 2010-2013 Philipp Benner
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
dnl we have a README.md instead of a README
AM_INIT_AUTOMAKE([gnu dist-bzip2 foreign])
dnl find c++ compiler but
dnl turn of AC_PROG_CC's default behaviour
: ${CXXFLAGS=""}
AC_PROG_CXX
dnl run tests (ac_try_compile) for C++ only
AC_LANG([C++])
AC_ISC_POSIX
AC_PROG_INSTALL
AM_PROG_CC_C_O
AC_PROG_LN_S
AM_PROG_LEX
AC_PROG_YACC
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL
LT_INIT
LT_PATH_LD

dnl ,---------------------------- 
dnl | GENERAL PACKAGE INFO
dnl `----------------------------

AC_SUBST([PACKAGE_SUMMARY], [""])
AC_SUBST([PACKAGE_URL], [""])
AC_SUBST([PACKAGE_LICENSE], ["GPL"])
AC_SUBST([PACKAGE_DESCRIPTION], [""])

dnl ,---------------------------- 
dnl | LIBRARIES
dnl `----------------------------

AC_CHECK_LIB([m],[cos],,
  AC_MSG_ERROR([this library is required]))
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,
  AC_MSG_ERROR([this library is required]))
AC_CHECK_LIB([gsl], [gsl_blas_dgemm],,
  AC_MSG_ERROR([this library is required]))
AC_CHECK_LIB([glpk], [glp_simplex], LIB_GLPK="-lglpk",
  AC_MSG_ERROR([this library is required]))
AC_CHECK_LIB([pthread], [pthread_create], LIB_PTHREAD="-lpthread",
  AC_MSG_ERROR([this library is required]))

AC_SUBST(LIB_GLPK)
AC_SUBST(LIB_PTHREAD)

dnl ,---------------------------- 
dnl | PYTHON
dnl `----------------------------

AM_PATH_PYTHON(2.5)
AC_PYTHON_DEVEL

dnl ,---------------------------- 
dnl | BOOST
dnl `----------------------------

AX_BOOST_BASE([1.47.0],,
  AC_MSG_ERROR([this library is required]))
AX_BOOST_REGEX(,
  AC_MSG_ERROR([this library is required]))
AX_BOOST_SYSTEM(,
  AC_MSG_ERROR([this library is required]))
AX_BOOST_THREAD(,
  AC_MSG_ERROR([this library is required]))
AX_BOOST_PYTHON(,
  AC_MSG_ERROR([this library is required]))

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
dnl | GCC_OPTIONS
dnl `----------------------------

AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
                [turn on debugging @<:@default=no@:>@])],
  [],
  [enable_debug=no])
AC_ARG_ENABLE([profiling],
  [AS_HELP_STRING([--enable-profiling],
                [turn on profiling support @<:@default=no@:>@])],
  [],
  [enable_profiling=no])
AC_ARG_ENABLE([profile-generate],
  [AS_HELP_STRING([--enable-profile-generate],
                [generate optimization profiles @<:@default=no@:>@])],
  [],
  [enable_profile_generate=no])
AC_ARG_ENABLE([lto],
  [AS_HELP_STRING([--enable-lto],
                [enable link time optimization (LTO) @<:@default=no@:>@])],
  [],
  [enable_lto=no])
AC_ARG_ENABLE([profile-use],
  [AS_HELP_STRING([--enable-profile-use],
                [use optimization profiles @<:@default=no@:>@])],
  [],
  [enable_profile_use=no])
AC_ARG_ENABLE([native],
  [AS_HELP_STRING([--enable-native],
                [compile with native support @<:@default=no@:>@])],
  [],
  [enable_native=no])
AS_IF([test "x$enable_debug" = "xyes"],
  [AX_CXXFLAGS_GCC_OPTION([-O0])
   AX_CXXFLAGS_GCC_OPTION([-g])
   AX_CXXFLAGS_GCC_OPTION([-ggdb])
   CXXFLAGS="$CXXFLAGS -DDEBUG"
   LDFLAGS="$LDFLAGS -static"],
  [AX_CXXFLAGS_GCC_OPTION([-O3])
   AX_CXXFLAGS_GCC_OPTION([-ffast-math])])
AS_IF([test "x$enable_lto" = "xyes"],
  [AX_CXXFLAGS_GCC_OPTION([-flto])])
AS_IF([test "x$enable_profiling" = "xyes"],
  [AX_CXXFLAGS_GCC_OPTION([-pg])
   LDFLAGS="$LDFLAGS -static"])
AS_IF([test "x$enable_profile_generate" = "xyes"],
  [AX_CXXFLAGS_GCC_OPTION([-fprofile-generate])
   LDFLAGS="$LDFLAGS -static"])
AS_IF([test "x$enable_profile_use" = "xyes"],
  [AX_CXXFLAGS_GCC_OPTION([-fprofile-use])
   LDFLAGS="$LDFLAGS"])
AS_IF([test "x$enable_native" = "xyes"],
  [AX_CXXFLAGS_GCC_OPTION([-march=native])])

AX_CXXFLAGS_GCC_OPTION([-std=c++11])
AX_CXXFLAGS_GCC_OPTION([-Wall])
AX_CXXFLAGS_GCC_OPTION([-Wwrite-strings])
AX_CXXFLAGS_GCC_OPTION([-Wstrict-prototypes])
AX_CXXFLAGS_GCC_OPTION([-Wno-trigraphs])
AX_CXXFLAGS_GCC_OPTION([-Wtrampolines])
#AX_CXXFLAGS_GCC_OPTION([-fPIC])
AX_CXXFLAGS_GCC_OPTION([-fno-nested-functions])
AX_CXXFLAGS_GCC_OPTION([-fomit-frame-pointer])
AX_CXXFLAGS_GCC_OPTION([-fvisibility-inlines-hidden])
dnl no runtime type information
dnl (-fno-rtti causes trouble while linking with rtti-enabled libraries)
#AX_CXXFLAGS_GCC_OPTION([-fno-rtti], CXXFLAGS_NORTTI)
dnl but it is possible to disable rtti in boost
AX_CXXFLAGS_GCC_OPTION([-DBOOST_NO_TYPEID], CXXFLAGS_NORTTI)
AX_CXXFLAGS_GCC_OPTION([-DBOOST_NO_RTTI], CXXFLAGS_NORTTI)
dnl turn link time optimization on (works only on some systems so far)
#AX_CXXFLAGS_GCC_OPTION([-flto])

# check for attributes
AX_GCC_ATTRIBUTE_HOT
AX_GCC_ATTRIBUTE_UNUSED

AC_SUBST(CXXFLAGS_NORTTI)

dnl ,---------------------------- 
dnl | AC_OUTPUT
dnl `----------------------------

AC_CONFIG_HEADER([tfbayes/config.h:tfbayes/config.h.in])
AC_CONFIG_FILES([
  Makefile
  data/Makefile
  tfbayes/version.h
  tfbayes/Makefile
  tfbayes/alignment/Makefile
  tfbayes/alignio/Makefile
  tfbayes/fastarithmetics/Makefile
  tfbayes/fasta/Makefile
  tfbayes/fg/Makefile
  tfbayes/config/Makefile
  tfbayes/interface/Makefile
  tfbayes/phylotree/Makefile
  tfbayes/dpm/Makefile
  tfbayes/cluster/Makefile
  tfbayes/uipac/Makefile
  tfbayes/weblogo/Makefile
  tfbayes/corebio/Makefile
  tfbayes/tools/Makefile
  tfbayes/utility/Makefile
  ])

AC_OUTPUT

dnl ,---------------------------- 
dnl | Print a summary
dnl `----------------------------

AC_MSG_NOTICE_BOX([Configure summary

System type:              ${host}
Build with debug symbols: ${enable_debug}
Profiling support:        ${enable_profiling}
Native processor support: ${enable_native}
])

dnl Local Variables:
dnl indent-tabs-mode: nil
dnl tab-stop-list: (2 4 6 8 10 12 14 16 18 20)
dnl End:
