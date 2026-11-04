#ifndef EQUATIONS_H
#define EQUATIONS_H

#include <complex.h>

/* Formula (2.2): intentionally direct, not numerically stable. */
void eq2(double a, double b, double c, double *x1, double *x2);

/* Formulas (2.3)-(2.4): stable calculation of real roots. */
void eq2_correct(double a, double b, double c, double *x1, double *x2);

/* Extension of eq2_correct for negative discriminants (C99 complex numbers). */
void eq2_correct_complex(double a, double b, double c,
                         double complex *x1, double complex *x2);

/* Formulas (2.6)-(2.13): x^3 + a*x^2 + b*x + c = 0. */
void eq3(double a, double b, double c,
         double complex *x1, double complex *x2, double complex *x3);

#endif
