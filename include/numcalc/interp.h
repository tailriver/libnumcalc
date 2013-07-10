#pragma once
#ifndef LIBNUMCALC_INTERP_H
#define LIBNUMCALC_INTERP_H

#include "typedef.h"

void interp1(
    integer nt,
    const double* xt,
    const double* yt,
    integer n,
    const double* x,
    double* y,
    integer* hint);

#endif /* LIBNUMCALC_INTERP_H */
