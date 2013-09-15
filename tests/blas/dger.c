#include <numcalc/blas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../common.h"

int main(int argc, char** argv) {
  double alpha, *x, *y, *A, *B, diff, difftemp;
  int m, n, incx, incy, lda, mm, nn, i, j, ix, iy;

  m = atoi(argv[1]);
  n = atoi(argv[2]);
  alpha = atof(argv[3]);
  incx = atoi(argv[4]);
  incy = atoi(argv[5]);

  mm = fabs(incx) * (m - 1) + 1;
  nn = fabs(incy) * (n - 1) + 1;
  lda = m + 1;

  x = (double*) malloc(sizeof(double) * mm);
  y = (double*) malloc(sizeof(double) * nn);
  A = (double*) malloc(sizeof(double) * lda * n);
  B = (double*) malloc(sizeof(double) * lda * n);

  for (i = 0; i < mm; i++) {
    x[i] = (double)(rand()) / RAND_MAX;
  }
  for (i = 0; i < nn; i++) {
    y[i] = (double)(rand()) / RAND_MAX;
  }
  for (i = 0; i < lda * n; i++) {
    A[i] = (double)(rand()) / RAND_MAX;
  }

  ix = incx > 0 ? 0 : (mm - 1);
  for (i = 0; i < m; i++, ix += incx) {
    iy = incy > 0 ? 0 : (nn - 1);
    for (j = 0; j < n; j++, iy += incy) {
      B[i+lda*j] = alpha * x[ix] * y[iy] + A[i+lda*j];
    }
  }

  dger(m, n, alpha, x, incx, y, incy, A, lda); 

  diff = 0.0;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      difftemp = fabs(A[i+lda*j] - B[i+lda*j]);
      if (difftemp > diff) diff = difftemp;
    }
  }
  printf("max diff. %e\n", diff);

  free(x);
  free(y);
  free(A);
  free(B);
}
