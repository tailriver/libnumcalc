#include "common.h"

double dasum(
    integer n,
    const double* dx,
    integer incx)
{
    integer i, ix;
    double asum;

    /* error for n < 0 */
    if (n < 1) return 0.0;

    /* fast return for unrollable loop */
    if (n == 1 || incx == 0) return *dx;

    asum = 0.0;
    if (incx == 1) {
        for (i = 0; i < n; i++) {
            asum += ABS(dx[i]);
        }
    } else {
        ix = incx > 0 ? 0 : -(n - 1) * incx;
        for (i = 0; i < n; i++) {
            asum += ABS(dx[ix]);
            ix += incx;
        }
    }
    return asum;
}
