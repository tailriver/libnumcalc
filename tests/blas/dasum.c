#include <numcalc/blas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../common.h"

int main(int argc, const char* argv[]) {
    double expected, got;
    double *x;
    int n, incx, nx, i;

    if (argc != 3) {
        fprintf(stderr, "Usage: ./dasum n incx\n");
        return EXIT_FAILURE;
    }

    n = atoi(argv[1]);
    incx = atoi(argv[2]);
    nx = (n - 1) * fabs(incx) + 1;

    x = malloc_double_array(nx, -1, 1);

    expected = 0.0;
    if (incx > 0) {
        for (i = 0; i < nx; i += incx) {
            expected += fabs(x[i]);
        }
    }
    else if (incx < 0) {
        for (i = nx - 1; i >= 0; i += incx) {
            expected += fabs(x[i]);
        }
    }
    else {
        expected = fabs(x[0]);
    }
    got = dasum(n, x, incx);
    printf("n=%d, incx=%d, expected: %e, got: %e\n", n, incx, expected, got);

    return EXIT_SUCCESS;
}
