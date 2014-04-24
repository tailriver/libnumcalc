#pragma once
#ifndef LIBNUMCALC_BLASEXT_H
#define LIBNUMCALC_BLASEXT_H

#include "typedef.h"

float ssum(
    integer n,
    const float* dx,
    integer incx);

double dsum(
    integer n,
    const double* dx,
    integer incx);

#endif /* LIBNUMCALC_BLASEXT_H */
