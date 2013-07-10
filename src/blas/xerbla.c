#include "common.h"

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif


void xerbla(
    const char* SRNAME,
    integer INFO)
{
#ifdef HAVE_STDIO_H
    fprintf(stderr, " ** On entry to %s parameter number %ld had an illegal value\n", SRNAME, INFO);
#endif
#ifdef HAVE_STDLIB_H
    exit(EXIT_FAILURE);
#endif
}
