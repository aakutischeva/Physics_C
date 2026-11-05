#ifndef EI_H
#define EI_H

#include <stddef.h>
#include <stdio.h>

/* Keep interactive table generation bounded to a reasonable file size. */
#define EI_MAX_POINTS 1000000U

/* Exponential integral Ei(x), defined for x > 0 in this assignment.
   Invalid arguments return NAN; Ei(+INFINITY) returns +INFINITY. */
double Ei(double x);

/* Write two whitespace-separated columns x and Ei(x) without a header.
   Count=1 writes only xmin. Returns 0 on success and -1 on error. */
int ei_write_table(FILE *stream, double xmin, double xmax, size_t count);

#endif
