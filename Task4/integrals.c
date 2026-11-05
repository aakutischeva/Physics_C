#include "integrals.h"

#include <errno.h>
#include <math.h>
#include <stddef.h>

typedef struct {
    long double sum;
    long double correction;
} Accumulator;

/* Neumaier-compensated summation for weakly cancelling terms. */
static void add(Accumulator *acc, long double x)
{
    long double next = acc->sum + x;
    if (fabsl(acc->sum) >= fabsl(x))
        acc->correction += (acc->sum - next) + x;
    else
        acc->correction += (x - next) + acc->sum;
    acc->sum = next;
}

static long double sum(const Accumulator *acc)
{
    return acc->sum + acc->correction;
}

static double invalid(int err)
{
    errno = err;
    return NAN;
}

static int valid_function_and_bounds(finteg f, double a, double b)
{
    return f != NULL && isfinite(a) && isfinite(b);
}

/* An integrand must return a finite double at every sampled node. */
static int evaluate(finteg f, long double x, long double *y)
{
    double argument = (double)x;
    if (!isfinite(argument)) {
        errno = ERANGE;
        return 0;
    }

    double value = f(argument);
    if (!isfinite(value)) {
        errno = EDOM;
        return 0;
    }

    *y = (long double)value;
    return 1;
}

/* Detect both overflow and underflow of a nonzero integral on return. */
static double checked_result(long double value)
{
    double converted = (double)value;
    if (!isfinite(value) || !isfinite(converted)
        || (value != 0.0L && converted == 0.0))
        return invalid(ERANGE);
    return converted;
}

/* Compute an abscissa using the wider type, with exact endpoints. */
static long double node(double a, double b, long double width,
                        unsigned long i, unsigned long n)
{
    if (i == 0UL)
        return (long double)a;
    if (i == n)
        return (long double)b;
    return (long double)a + width * ((long double)i / (long double)n);
}

double trapz(finteg f, double a, double b, int n)
{
    if (!valid_function_and_bounds(f, a, b) || n <= 0)
        return invalid(EINVAL);
    errno = 0;

    if (a == b)
        return 0.0;

    long double width = (long double)b - (long double)a;
    long double h = width / (long double)n;
    if (!isfinite(width) || !isfinite(h) || h == 0.0L)
        return invalid(ERANGE);

    Accumulator acc = {0.0L, 0.0L};
    unsigned long intervals = (unsigned long)n;
    for (unsigned long i = 0UL; i <= intervals; ++i) {
        long double y;
        if (!evaluate(f, node(a, b, width, i, intervals), &y))
            return NAN;
        add(&acc, (i == 0UL || i == intervals) ? 0.5L * y : y);
    }

    return checked_result(h * sum(&acc));
}

double simpson(finteg f, double a, double b, int n)
{
    if (!valid_function_and_bounds(f, a, b) || n < 5)
        return invalid(EINVAL);
    errno = 0;

    if (a == b)
        return 0.0;

    long double width = (long double)b - (long double)a;
    long double h = width / (long double)n;
    if (!isfinite(width) || !isfinite(h) || h == 0.0L)
        return invalid(ERANGE);

    Accumulator acc = {0.0L, 0.0L};
    unsigned long intervals = (unsigned long)n;
    for (unsigned long i = 0UL; i <= intervals; ++i) {
        long double y;
        if (!evaluate(f, node(a, b, width, i, intervals), &y))
            return NAN;

        long double weight;
        if (i == 0UL || i == intervals)
            weight = 3.0L / 8.0L;
        else if (i == 1UL || i == intervals - 1UL)
            weight = 7.0L / 6.0L;
        else if (i == 2UL || i == intervals - 2UL)
            weight = 23.0L / 24.0L;
        else
            weight = 1.0L;

        add(&acc, weight * y);
    }

    return checked_result(h * sum(&acc));
}

double trapz_adapt(finteg f, double a, double b, double eps)
{
    if (!valid_function_and_bounds(f, a, b) || !isfinite(eps)
        || eps <= 0.0 || eps >= 1.0)
        return invalid(EINVAL);
    errno = 0;

    if (a == b)
        return 0.0;

    long double width = (long double)b - (long double)a;
    if (!isfinite(width))
        return invalid(ERANGE);

    long double fa, fb;
    if (!evaluate(f, (long double)a, &fa)
        || !evaluate(f, (long double)b, &fb))
        return NAN;

    long double previous = width * (fa + fb) / 2.0L;
    /* The stopping scale is the trapezoidal approximation to
       integral_a^b |f(x)| |dx|.  It stays positive when bounds are
       reversed and does not disappear through sign cancellation. */
    long double absolute_previous =
        fabsl(width) * (fabsl(fa) + fabsl(fb)) / 2.0L;
    if (!isfinite(previous) || !isfinite(absolute_previous))
        return invalid(ERANGE);

    unsigned long intervals = 1UL;
    /* 22 refinements mean at most 2^22 new midpoint evaluations. */
    for (int iteration = 0; iteration < 22; ++iteration) {
        long double h = width / (long double)intervals;
        if (!isfinite(h) || h == 0.0L)
            return invalid(ERANGE);

        Accumulator midpoints = {0.0L, 0.0L};
        Accumulator absolute_midpoints = {0.0L, 0.0L};
        for (unsigned long i = 0UL; i < intervals; ++i) {
            long double x = (long double)a + width *
                ((long double)i + 0.5L) / (long double)intervals;
            long double y;
            if (!evaluate(f, x, &y))
                return NAN;
            add(&midpoints, y);
            add(&absolute_midpoints, fabsl(y));
        }

        /* T(2n) = T(n)/2 + h/2 * sum_i f(a + (i+1/2) h). */
        long double current = previous / 2.0L + h * sum(&midpoints) / 2.0L;
        long double absolute_current = absolute_previous / 2.0L
            + fabsl(h) * sum(&absolute_midpoints) / 2.0L;
        if (!isfinite(current) || !isfinite(absolute_current))
            return invalid(ERANGE);

        /* Scale-relative stopping test.  A fixed max(1, |current|)
           would be much too loose for integrals with magnitude < 1. */
        long double tolerance = (long double)eps * absolute_current;
        if (fabsl(current - previous) <= tolerance)
            return checked_result(current);

        previous = current;
        absolute_previous = absolute_current;
        intervals *= 2UL;
    }

    return invalid(ERANGE); /* The refinement limit was reached. */
}
