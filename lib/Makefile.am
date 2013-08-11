AM_CPPFLAGS = -I$(top_srcdir)/include
AM_CFLAGS = -Wall
AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libnumcalc.a

libnumcalc_a_SOURCES = \
    blas/dasum.c \
    blas/daxpy.c \
    blas/dcopy.c \
    blas/dgemv.c \
    blas/dger.c \
    blas/dscal.c \
    blas/dswap.c \
    blas/dtbsv.c \
    blas/idamax.c \
    blas/lsame.c \
    blas/xerbla.c \
    blasext/dsum.c \
    interp/interp1.c \
    lapack/dgbtf2.c \
    lapack/dgbtrs.c

noinst_HEADERS = \
    array_notation.h \
    common.h \
    blas/common.h \
    blasext/common.h \
    interp/common.h \
    lapack/common.h
