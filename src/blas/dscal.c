#include "blas_private.h"


void dscal(
    integer n,
    double da,
    double* dx,
    integer incx)
{
    integer i, ix;

    /* fast return for n */
    if (n < 1) return;

    /* fast return for da */
    if (da == 1.0) return;

    /* fast return for unrollable loop */
    if (n == 1) {
        dx[0] *= da;
        return;
    }

    /* no fast return for incx = 0 */

    if (incx == 1) {
        for (i = 0; i < n; i++) {
            dx[i] *= da;
        }
    } else {
        ix = incx > 0 ? 0 : -(n - 1) * incx;
        for (i = 0; i < n; i++) {
            dx[ix] *= da;
            ix += incx;
        }
    }
}
