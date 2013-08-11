#include <numcalc/blas.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    double x[10], expected, got;
    size_t i;

    srand(time(NULL));
    for (i = 0; i < 10; i++) {
        x[i] = (double) (rand() - RAND_MAX / 2) / RAND_MAX;
    }

    /* incx = 1 */
    expected = 0.0;
    for (i = 0; i < 10; i++) {
        expected += fabs(x[i]);
    }
    got = dasum(10, x, 1);
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = 2 */
    expected = 0.0;
    for (i = 0; i < 10; i += 2) {
        expected += fabs(x[i]);
    }
    got = dasum(5, x, 2);
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = -1 */
    expected = 0.0;
    for (i = 9; i > 1; i--) {
        expected += fabs(x[i]);
    }
    got = dasum(8, &x[2], -1);
    printf("expected: %e, got: %e\n", expected, got);

    /* incx = -2 */
    expected = 0.0;
    for (i = 9; i > 2; i -= 2) {
        expected += fabs(x[i]);
    }
    got = dasum(4, &x[3], -2);
    printf("expected: %e, got: %e\n", expected, got);

    return EXIT_SUCCESS;
}
