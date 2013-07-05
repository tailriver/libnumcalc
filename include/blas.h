#pragma once
#include "typedef.h"

double dasum(
    integer n,
    double* dx,
    integer incx);


void daxpy(
    integer n,
    double da,
    const double* dx,
    integer incx,
    double* dy,
    integer incy);


void dcopy(
      integer n,
      const double* dx,
      integer incx,
      double* dy,
      integer incy);


void dgemv(
      char TRANS,
      integer M,
      integer N,
      double ALPHA,
      const double* A,
      integer LDA,
      const double* X,
      integer INCX,
      double BETA,
      double* Y,
      integer INCY);


void dger(
      integer M,
      integer N,
      double ALPHA,
      const double* X,
      integer INCX,
      const double* Y,
      integer INCY,
      double* A,
      integer LDA);


void dscal(
      integer N,
      double DA,
      double* DX,
      integer INCX);


void dswap(
      integer N,
      double* DX,
      integer INCX,
      double* DY,
      integer INCY);


void dtbsv(
      char UPLO,
      char TRANS,
      char DIAG,
      integer N,
      integer K,
      const double* A,
      integer LDA,
      double* X,
      integer INCX);


void dtrsv(
      char UPLO,
      char TRANS,
      char DIAG,
      integer N,
      const double* A,
      integer LDA,
      double* X,
      integer INCX);


integer idamax(
      integer N,
      const double* DX,
      integer INCX);


logical lsame(
      char CA,
      char CB);


void xerbla(
      const char* SRNAME,
      integer N);