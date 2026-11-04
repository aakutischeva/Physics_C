#ifndef EQUATIONS_H
#define EQUATIONS_H

#include <complex.h>

void eq2(double a, double b, double c, double *x1, double *x2);

void eq2_correct(double a, double b, double c, double *x1, double *x2);

void eq2_correct_complex(double a, double b, double c,
                         double complex *x1, double complex *x2);

void eq3(double a, double b, double c,
         double complex *x1, double complex *x2, double complex *x3);

#endif
