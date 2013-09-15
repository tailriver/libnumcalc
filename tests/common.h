#include <stdlib.h>
#include <time.h>

static int srand_called = 0;

void check_srand_called() {
    if (srand_called) return;

    srandom(time(NULL));
    srand_called = 1;
}

double* malloc_double_array(int n, double min, double max) {
    double* x;
    int i;

    check_srand_called();

    x = (double*) malloc(sizeof(double) * n);
    for (i = 0; i < n; i++) {
        x[i] = min + ((double)random() / 2147483647) * (max - min);
    }

    return x;
}
