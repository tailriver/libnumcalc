#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <blas.h>

#define ABS(a) ((a) < 0 ? (-(a)) : (a))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define X(i) X[(i)-1]
#define Y(i) Y[(i)-1]

#define A(i,j) A[((i)-1)+((j)-1)*LDA]
#define B(i,j) B[((i)-1)+((j)-1)*LDB]
#define C(i,j) C[((i)-1)+((j)-1)*LDC]
#define AB(i,j) AB[((i)-1)+((j)-1)*LDAB]
#define PAB(i,j) &AB[((i)-1)+((j)-1)*LDAB]
