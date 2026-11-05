#ifndef LINAL_H
#define LINAL_H

/* Functions specified in section 2.5 of Marenkov's textbook. */
double **alloc2d(int n1, int n2);
void free2d(double **a, int n1, int n2);
double vabs(double *a, int n);
double scal_prod(double *a, double *b, int n);
void fill(double **a, int n1, int n2, double x);

/* The original API has only two dimensions. Its explicit contract here:
   a[n1][n2] * b[n2][n1] -> c[n1][n1]. */
void matrix_prod(double **a, double **b, double **c, int n1, int n2);
void show2d(double **a, int n1, int n2);

/* Extension for a general rectangular product:
   a[rows][inner] * b[inner][cols] -> c[rows][cols]. */
void matrix_prod_rect(double **a, double **b, double **c,
                      int rows, int inner, int cols);

/* Error convention: scalar functions return NAN; alloc2d returns NULL;
   void routines leave the output unchanged on invalid input whenever
   possible. Consult errno: EINVAL, EDOM, ERANGE, ENOMEM. */

#endif
