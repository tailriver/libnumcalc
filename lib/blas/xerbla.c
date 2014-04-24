#include "common.h"

#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif


void xerbla(
    const char* SRNAME,
    integer INFO)
{
    fprintf(stderr, " ** On entry to %s parameter number %ld had an illegal value\n", SRNAME, INFO);
#ifdef HAVE_STDLIB_H
    exit(EXIT_FAILURE);
#endif
}
