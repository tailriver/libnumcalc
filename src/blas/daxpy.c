#include "common.h"


void daxpy(
    integer n,
    double da,
    const double* dx,
    integer incx,
    double* dy,
    integer incy)
{
    integer i, ix, iy;

    /* fast return for n */
    if (n < 1) return;

    /* fast return for da */
    if (da == 0.0) return;

    if (incx == 0) {
        da *= *dx;
        iy = incy > 0 ? 0 : -(n - 1) * incy;
        for (i = 0; i < n; i++) {
            dy[iy] += da;
            iy += incy;
        }
    } else if (incx == 1 && incy == 1) {
        for (i = 0; i < n; i++) {
            dy[i] += da * dx[i];
        }
    } else {
        ix = incx > 0 ? 0 : -(n - 1) * incx;
        iy = incy > 0 ? 0 : -(n - 1) * incy;
        for (i = 0; i < n; i++) {
            dy[iy] += da * dx[ix];
            ix += incx;
            iy += incy;
        }
    }
}
