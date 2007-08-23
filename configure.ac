m4_define([ngila_rev], [m4_esyscmd([svnversion -n .])])
AC_INIT([ngila], [1.2.1-RELEASE], [http://scit.us/projects/ngila/])
AC_CONFIG_SRCDIR([src/ngila.h])
AC_CONFIG_AUX_DIR([config])
AC_PREREQ(2.56)
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([config.h])

dnl Program Options
AC_ARG_ENABLE(all-static,
	[  --enable-all-static     Force static linking],
	USE_ALL_STATIC="$enableval", USE_ALL_STATIC="no")

AC_ARG_ENABLE(debug,
	[  --enable-debug          Enable debugging information],
	USE_DEBUG="$enableval", USE_DEBUG="no")

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP

AC_USE_SYSTEM_EXTENSIONS

if test $USE_DEBUG = yes ; then 
  DEBUG=1
  CFLAGS="-Wall -g -O0"
  CXXFLAGS="-Wall -g -O0"
else
  DEBUG=0
  AX_CC_MAXOPT
  CFLAGS="$CFLAGS -DNDEBUG"
  CXXFLAGS="$CFLAGS"
fi

AC_PROG_LIBTOOL

if test $USE_ALL_STATIC = yes ; then
  ALL_STATIC_LDFLAGS="-all-static"
else
  ALL_STATIC_LDFLAGS=""
fi
AC_SUBST(ALL_STATIC_LDFLAGS)

AC_HEADER_STDC
AM_PATH_GSL([1.0], [], AC_MSG_ERROR(could not find required version of GSL))

AC_C_INLINE
if test "$ac_cv_c_inline" != no ; then
	AC_DEFINE(HAVE_INLINE, 1, [HAVE_INLINE])
	AC_SUBST(HAVE_INLINE)
fi

AC_LANG_SAVE
AC_LANG_CPLUSPLUS

AC_DEFUN([AC_EXT_CHECK_LIB],
[
	AX_EXT_HAVE_LIB([/lib /usr/lib /usr/local/lib /opt], [$1], [$2], [$5])
	AS_IF([test AS_VAR_GET([ext_lib_found]) = yes], [$3], [$4])
])

AX_BOOST_BASE([1.33.0])
AX_BOOST_PROGRAM_OPTIONS

dnl Checks for header files.
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
dnl AC_C_BIGENDIAN

dnl Checks for library functions.

AX_CXXFLAGS_WARN_ALL
AX_CFLAGS_WARN_ALL

AX_CREATE_STDINT_H

AC_LANG_RESTORE

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile doc/Makefile matrix/Makefile vs/Makefile])
AC_OUTPUT
