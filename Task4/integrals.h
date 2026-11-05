#ifndef INTEGRALS_H
#define INTEGRALS_H

/* Function-pointer type required by section 2.4 of the textbook. */
typedef double (*finteg)(double);

/* Formula (2.21): composite trapezoidal rule, n >= 1. */
double trapz(finteg f, double a, double b, int n);

/* Formula (2.23): fourth-order endpoint-corrected rule, n >= 5.
   The name simpson is prescribed by the textbook; this is NOT
   the alternating-weight composite Simpson 1/3 rule (2.22). */
double simpson(finteg f, double a, double b, int n);

/* Adaptive trapezoidal rule, 0 < eps < 1, up to 22 refinements.
   All routines return NAN on errors and set errno to EINVAL, EDOM,
   or ERANGE; finite results are returned as double. */
double trapz_adapt(finteg f, double a, double b, double eps);

#endif
