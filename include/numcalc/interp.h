#pragma once
#ifndef LIBNUMCALC_INTERP_H
#define LIBNUMCALC_INTERP_H

#include "typedef.h"

void sinterp1(
    integer nt,
    const float* xt,
    const float* yt,
    integer n,
    const float* x,
    float* y,
    integer* hint);

void dinterp1(
    integer nt,
    const double* xt,
    const double* yt,
    integer n,
    const double* x,
    double* y,
    integer* hint);

#endif /* LIBNUMCALC_INTERP_H */
