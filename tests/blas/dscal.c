#include <numcalc/blas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../common.h"

int main(int argc, const char* argv[]) {
    double *x, *y, da, diff;
    int n, incx, nx, i;

    if (argc != 4) {
        fprintf(stderr, "Usage: ./dscal n da incx\n");
        return EXIT_FAILURE;
    }

    n = atoi(argv[1]);
    da = atof(argv[2]);
    incx = atoi(argv[3]);
    nx = (n - 1) * fabs(incx) + 1;

    x = malloc_double_array(nx, -1, 1);
    y = malloc_double_array(nx, -1, 1);
    dcopy(nx, x, 1, y, 1);

    dscal(nx, da, x, incx);
    if (incx > 0) {
        for (i = 0; i < nx; i += incx) {
            diff += x[i] - da * y[i];
        }
    }
    else if (incx < 0) {
        for (i = nx - 1; i >= 0; i += incx) {
            diff += x[i] - da * y[i];
        }
    }
    else {
        i = 0;
        diff = x[i] - da * n * y[i];
    }
    printf("n=%d, da=%e, incx=%d, diff.: %e\n", n, da, incx, diff);

    return EXIT_SUCCESS;
}
