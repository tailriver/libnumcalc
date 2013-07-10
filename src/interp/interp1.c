#include "common.h"

void interp1(
    integer nt,
    const double* xt,
    const double* yt,
    integer n,
    const double* x,
    double* y,
    integer* hint)
{
    integer i, t;

    if (nt < 2 || n < 1) return;

    /* hint has left index (for user convenience) */
    /* t has right index */

    if (xt[0] < xt[1]) {

        /* xt is ascending order */
        for (i = 0; i < n; i++) {
            if (xt[nt-1] < x[i]) {
                t = nt - 1;
            } else {
                t = hint == NULL ? 1 : hint[i] + 1;
                for ( ; xt[t] < x[i] && t < nt; t++); /* empty loop */
            }
            y[i] = yt[t-1] + (x[i] - xt[t-1]) / (xt[t] - xt[t-1]) * (yt[t] - yt[t-1]);

            if (hint != NULL) hint[i] = t - 1;
        }

    } else {

        /* xt is descending order */
        for (i = 0; i < n; i++) {
            if (xt[nt-1] > x[i]) {
                t = nt - 1;
            } else {
                t = hint == NULL ? 1 : hint[i] + 1;
                for ( ; xt[t] > x[i] && t < nt; t++); /* empty loop */
            }
            y[i] = yt[t-1] + (x[i] - xt[t-1]) / (xt[t] - xt[t-1]) * (yt[t] - yt[t-1]);

            if (hint != NULL) hint[i] = t - 1;
        }
    }
}
