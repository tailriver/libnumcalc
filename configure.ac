#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([numcalc], [0.1.1], [shinsuke@tailriver.net])
AC_CONFIG_SRCDIR([lib/blas/daxpy.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([m], [sin])

# Checks for header files.
AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([time.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile
    include/Makefile
    include/numcalc/Makefile
    lib/Makefile
    tests/Makefile
    tests/blas/Makefile
    tests/blasext/Makefile
    tests/interp/Makefile
])
AC_OUTPUT
