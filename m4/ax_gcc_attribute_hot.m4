dnl Copyright (C) 2013 Philipp Benner
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

AC_DEFUN([AX_GCC_ATTRIBUTE_HOT],[
AC_CACHE_CHECK(
  [for  __attribute__((hot))],
  ax_cv_gcc_attribute_hot,[
  ac_save_cxx_werror_flag=$ac_cxx_werror_flag
  ac_cxx_werror_flag=yes
  AC_TRY_COMPILE([
    __attribute__((hot))
    int f(int i) { return i; }],
    [],
    ax_cv_gcc_attribute_hot=yes, ax_cv_gcc_attribute_hot=no)
  ac_cxx_werror_flag=$ac_save_cxx_werror_flag
  ])
  if test "$ax_cv_gcc_attribute_hot" = yes; then
    AC_DEFINE([GCC_ATTRIBUTE_HOT],[__attribute__((hot))],
      [Attribute to declare functions as hot spots of the program])
  else
    AC_DEFINE([GCC_ATTRIBUTE_HOT],[],
      [Attribute to declare functions as hot spots of the program])
  fi
])
