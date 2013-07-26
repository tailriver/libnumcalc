#include "common.h"

void dswap(
    integer n,
    double* dx,
    integer incx,
    double* dy,
    integer incy)
{
    integer i, ix, iy;
    double temp;

    /* fast return for n */
    if (n < 1) return;

    /* fast return for same pointer */
    if (dx == dy && incx == incy) return;

    if (incx == 1 && incy == 1) {
        for (i = 0; i < n; i++) {
            temp = dx[i];
            dx[i] = dy[i];
            dy[i] = temp;
        }
    } else {
        ix = incx > 0 ? 0 : -(n - 1) * incx;
        iy = incy > 0 ? 0 : -(n - 1) * incy;
        for (i = 0; i < n; i++) {
            temp = dx[ix];
            dx[ix] = dy[iy];
            dy[iy] = temp;
            ix += incx;
            iy += incy;
        }
    }
}
