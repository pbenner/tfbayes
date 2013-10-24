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

AC_DEFUN([AX_BOOST_REGEX],
[
AC_CHECK_LIB(boost_regex-gcc-mt, main, BOOST_REGEX_LIB="-lboost_regex-gcc-mt",
  [AC_CHECK_LIB(boost_regex-mt, main, BOOST_REGEX_LIB="-lboost_regex-mt",
    [AC_CHECK_LIB(boost_regex, main, BOOST_REGEX_LIB="-lboost_regex",
      [$2])])])

AC_SUBST(BOOST_REGEX_LIB)
])
