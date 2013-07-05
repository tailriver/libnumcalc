#include "blas_private.h"

#ifdef HAVE_STRING_H
#include <string.h>
#endif


void dcopy(
    integer n,
    const double* dx,
    integer incx,
    double* dy,
    integer incy)
{
    integer i, ix, iy;

    /* fast return for n */
    if (n < 1) return;

    /* fast return for same pointer */
    if (dx == dy && incx == incy) return;

    if (incx == 0 && incy == 1) {
        for (i = 0; i < n; i++) {
            dy[i] = *dx;
        }
    } else if (incx == 1 && incy == 1) {
#ifdef HAVE_STRING_H
        memcpy(dy, dx, sizeof(double) * n);
#else
        for (i = 0; i < n; i++) {
            dy[i] = dx[i];
        }
#endif
    } else {
        ix = incx > 0 ? 0 : -(n - 1) * incx;
        iy = incy > 0 ? 0 : -(n - 1) * incy;
        for (i = 0; i < n; i++) {
            dy[iy] = dx[ix];
            ix += incx;
            iy += incy;
        }
    }
}
