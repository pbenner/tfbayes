dnl -- acsite.m4 --
dnl
dnl Copyright (C) 2002, 2003, 2004, 2006 Philipp Benner
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

dnl ,------------------------------------
dnl | MACRO DEFINITIONS
dnl `------------------------------------

AC_DEFUN([AC_MSG_WARN_BOX],
  [# display a nice ascii message box
	  AC_MSG_WARN([,---------------])
	  echo -e "$1" | while read LINE;
		  do AC_MSG_WARN([| $LINE]);
		done
	  AC_MSG_WARN([`---------------])])

AC_DEFUN([AC_MSG_NOTICE_BOX],
  [# display a nice ascii message box
	  AC_MSG_NOTICE([,---------------])
	  echo -e "$1" | while read LINE;
		  do AC_MSG_NOTICE([| $LINE]);
		done
	  AC_MSG_NOTICE([`---------------])])

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

dnl Local Variables:
dnl tab-width: 2
dnl indent-tabs-mode: t
dnl End:
dnl vim: tabstop=2 noexpandtab
