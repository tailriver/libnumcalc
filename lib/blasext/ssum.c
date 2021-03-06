#include "common.h"

float ssum(
    integer n,
    const float* dx,
    integer incx)
{
    integer i, ix;
    float sum;

    /* error for n < 0 */
    if (n < 1) return 0.0;

    /* fast return for unrollable loop */
    if (n == 1 || incx == 0) return *dx;

    sum = 0.0;
    if (incx == 1) {
        for (i = 0; i < n; i++) {
            sum += dx[i];
        }
    } else {
        ix = incx > 0 ? 0 : -(n - 1) * incx;
        for (i = 0; i < n; i++) {
            sum += dx[ix];
            ix += incx;
        }
    }
    return sum;
}
