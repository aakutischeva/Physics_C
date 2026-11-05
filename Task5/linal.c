#include "linal.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    long double sum;
    long double correction;
} Accumulator;

static void add(Accumulator *acc, long double value)
{
    long double next = acc->sum + value;
    if (fabsl(acc->sum) >= fabsl(value))
        acc->correction += (acc->sum - next) + value;
    else
        acc->correction += (value - next) + acc->sum;
    acc->sum = next;
}

static long double total(const Accumulator *acc)
{
    return acc->sum + acc->correction;
}

static double checked_double(long double value)
{
    double result = (double)value;
    if (!isfinite(value) || !isfinite(result)
        || (value != 0.0L && result == 0.0)) {
        errno = ERANGE;
        return NAN;
    }
    return result;
}

static int valid_rows(double **a, int rows)
{
    if (a == NULL || rows <= 0)
        return 0;
    for (int i = 0; i < rows; ++i)
        if (a[i] == NULL)
            return 0;
    return 1;
}

double **alloc2d(int n1, int n2)
{
    if (n1 <= 0 || n2 <= 0) {
        errno = EINVAL;
        return NULL;
    }

    size_t rows = (size_t)n1;
    size_t cols = (size_t)n2;
    if (rows > SIZE_MAX / sizeof(double *)
        || cols > SIZE_MAX / sizeof(double)) {
        errno = ERANGE;
        return NULL;
    }

    double **matrix = calloc(rows, sizeof(*matrix));
    if (matrix == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    for (size_t i = 0; i < rows; ++i) {
        matrix[i] = calloc(cols, sizeof(**matrix));
        if (matrix[i] == NULL) {
            for (size_t k = 0; k < i; ++k)
                free(matrix[k]);
            free(matrix);
            errno = ENOMEM;
            return NULL;
        }
    }

    errno = 0;
    return matrix;
}

void free2d(double **a, int n1, int n2)
{
    if (a == NULL)
        return;
    if (n1 <= 0 || n2 <= 0) {
        errno = EINVAL;
        return;
    }

    for (int i = 0; i < n1; ++i)
        free(a[i]);
    free(a);
    errno = 0;
}

/* Scaled sum of squares: avoids squaring a huge or tiny component. */
double vabs(double *a, int n)
{
    if (a == NULL || n <= 0) {
        errno = EINVAL;
        return NAN;
    }

    long double scale = 0.0L;
    long double ssq = 1.0L;
    for (int i = 0; i < n; ++i) {
        if (!isfinite(a[i])) {
            errno = EDOM;
            return NAN;
        }
        long double magnitude = fabsl((long double)a[i]);
        if (magnitude != 0.0L) {
            if (scale < magnitude) {
                long double ratio = scale / magnitude;
                ssq = 1.0L + ssq * ratio * ratio;
                scale = magnitude;
            } else {
                long double ratio = magnitude / scale;
                ssq += ratio * ratio;
            }
        }
    }

    errno = 0;
    return checked_double(scale == 0.0L ? 0.0L : scale * sqrtl(ssq));
}

double scal_prod(double *a, double *b, int n)
{
    if (a == NULL || b == NULL || n <= 0) {
        errno = EINVAL;
        return NAN;
    }

    Accumulator acc = {0.0L, 0.0L};
    for (int i = 0; i < n; ++i) {
        if (!isfinite(a[i]) || !isfinite(b[i])) {
            errno = EDOM;
            return NAN;
        }
        add(&acc, (long double)a[i] * (long double)b[i]);
    }

    errno = 0;
    return checked_double(total(&acc));
}

void fill(double **a, int n1, int n2, double x)
{
    if (n2 <= 0 || !valid_rows(a, n1)) {
        errno = EINVAL;
        return;
    }
    if (!isfinite(x)) {
        errno = EDOM;
        return;
    }

    for (int i = 0; i < n1; ++i)
        for (int j = 0; j < n2; ++j)
            a[i][j] = x;
    errno = 0;
}

void matrix_prod_rect(double **a, double **b, double **c,
                      int rows, int inner, int cols)
{
    if (rows <= 0 || inner <= 0 || cols <= 0
        || !valid_rows(a, rows) || !valid_rows(b, inner)
        || !valid_rows(c, rows)) {
        errno = EINVAL;
        return;
    }

    /* Compute into a temporary result first.  This ensures correct
       results for c == a/b when the shapes permit it, and prevents
       partial writes after a numeric error. */
    double **work = alloc2d(rows, cols);
    if (work == NULL)
        return;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Accumulator acc = {0.0L, 0.0L};
            for (int k = 0; k < inner; ++k) {
                if (!isfinite(a[i][k]) || !isfinite(b[k][j])) {
                    free2d(work, rows, cols);
                    errno = EDOM;
                    return;
                }
                add(&acc, (long double)a[i][k] * (long double)b[k][j]);
            }
            work[i][j] = checked_double(total(&acc));
            if (!isfinite(work[i][j])) {
                free2d(work, rows, cols);
                errno = ERANGE;
                return;
            }
        }
    }

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            c[i][j] = work[i][j];
    free2d(work, rows, cols);
    errno = 0;
}

void matrix_prod(double **a, double **b, double **c, int n1, int n2)
{
    matrix_prod_rect(a, b, c, n1, n2, n1);
}

void show2d(double **a, int n1, int n2)
{
    if (n2 <= 0 || !valid_rows(a, n1)) {
        errno = EINVAL;
        return;
    }

    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            if (printf("% .9g%s", a[i][j], j + 1 == n2 ? "\n" : " ") < 0) {
                if (errno == 0)
                    errno = EIO;
                return;
            }
        }
    }
    errno = 0;
}
