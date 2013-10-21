AC_DEFUN([AX_GCC_ATTRIBUTE_HOT],[dnl
AC_CACHE_CHECK(
 [the compiler for  __attribute__((hot))],
 ax_cv_gcc_const_call,[
 AC_TRY_COMPILE([__attribute__((hot))
 int f(int i) { return i; }],
 [],
 ax_cv_gcc_attribute_hot=yes, ax_cv_gcc_attribute_hot=no)])
 if test "$ax_cv_gcc_attribute_hot" = yes; then
   AC_DEFINE([GCC_ATTRIBUTE_HOT],[__attribute__((hot))],
    [most gcc compilers know a function __attribute__((hot))])
 fi
])
