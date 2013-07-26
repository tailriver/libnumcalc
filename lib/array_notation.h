#pragma once

/* vector */
#define X(i) X[(i)-1]
#define Y(i) Y[(i)-1]

/* matrix */
#define A(i,j) A[((i)-1)+((j)-1)*LDA]
#define B(i,j) B[((i)-1)+((j)-1)*LDB]
#define C(i,j) C[((i)-1)+((j)-1)*LDC]

/* banded matrix */
#define AB(i,j) AB[((i)-1)+((j)-1)*LDAB]

/* pivot matrix */
#define IPIV(i) IPIV[(i)-1]
