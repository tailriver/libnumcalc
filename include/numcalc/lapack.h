#pragma once
#ifndef LIBNUMCALC_LAPACK_H
#define LIBNUMCALC_LAPACK_H

#include "typedef.h"

void dgbtf2(
      integer M,
      integer N,
      integer KL,
      integer KU,
      double* AB,
      integer LDAB,
      integer* IPIV,
      integer* INFO);

void dgbtrs(
      char TRANS,
      integer N,
      integer KL,
      integer KU,
      integer NRHS,
      const double* AB,
      integer LDAB,
      const integer* IPIV,
      double* B,
      integer LDB,
      integer* INFO);

#endif /* LIBNUMCALC_LAPACK_H */
