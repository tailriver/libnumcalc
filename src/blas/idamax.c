#include "blas_private.h"

/* returns zero-based index */

integer idamax(
    integer n,
    const double* dx,
    integer incx)
{
    integer i, ix;
    integer iamax;
    double amax, temp;

    /* error for n < 0 */
    if (n < 1) return -1;

    /* fast return for unrollable loop */
    if (n == 1 || incx == 0) return 0;

    if (incx == 1) {
        iamax = 0;
        amax = ABS(dx[0]);
        for (i = 1; i < n; i++) {
            temp = ABS(dx[i]);
            if (temp > amax) {
                amax = temp;
                iamax = i;
            }
        }
    } else {
        ix = incx > 0 ? 0 : -(n - 1) * incx;
        iamax = ix;
        amax = ABS(dx[ix]);
        ix += incx;
        for (i = 1; i < n; i++) {
            temp = ABS(dx[ix]);
            if (temp > amax) {
                amax = temp;
                iamax = i;
            }
            ix += incx;
        }
    }
    return iamax;
}
