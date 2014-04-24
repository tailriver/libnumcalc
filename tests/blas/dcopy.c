#include <numcalc/blas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void clear(size_t n, double* x) {
    size_t i;

    for (i = 0; i < n; i++) {
        x[i] = sqrt(-1);
    }
}

int main(void) {
    double x[10], y[10], expected, got;
    size_t i;

    srand(time(NULL));
    for (i = 0; i < 10; i++) {
        x[i] = (double) (rand() - RAND_MAX / 2) / RAND_MAX;
    }
    expected = 0.0;

    /* incx = 1, incy = 1 */
    clear(10, y);
    dcopy(10, x, 1, y, 1);
    got = 0.0;
    for (i = 0; i < 10; i++) {
        got += x[i] - y[i];
    }
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = 2, incy = 1 */
    clear(10, y);
    dcopy(5, x, 2, y, 1);
    got = 0.0;
    for (i = 0; i < 5; i++) {
        got += x[2*i] - y[i];
    }
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = 1, incy = 2 */
    clear(10, y);
    dcopy(5, x, 1, y, 2);
    got = 0.0;
    for (i = 0; i < 5; i++) {
        got += x[i] - y[2*i];
    }
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = 2, incy = 2 */
    clear(10, y);
    dcopy(5, x, 2, y, 2);
    got = 0.0;
    for (i = 0; i < 5; i++) {
        got += x[2*i] - y[2*i];
    }
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = 1, incy = -1 */
    clear(10, y);
    dcopy(10, x, 1, y, -1);
    got = 0.0;
    for (i = 0; i < 10; i++) {
        got += x[i] - y[9-i];
    }
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = 1, incy = -2 */
    clear(10, y);
    dcopy(5, x, 1, y, -2);
    got = 0.0;
    for (i = 0; i < 5; i++) {
        got += x[i] - y[8-2*i];
    }
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = -1, incy = 1 */
    clear(10, y);
    dcopy(10, x, -1, y, 1);
    got = 0.0;
    for (i = 0; i < 10; i++) {
        got += x[9-i] - y[i];
    }
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = -2, incy = 1 */
    clear(10, y);
    dcopy(5, x, -2, y, 1);
    got = 0.0;
    for (i = 0; i < 5; i++) {
        got += x[8-2*i] - y[i];
    }
    printf("expected: %e, got: %e\n", expected, got);

    return EXIT_SUCCESS;
}
