#include "ei.h"

#include <math.h>

#define EULER 0.57721566490153286061L
#define MAXIT 100
#define EPS 1e-7L

/* (2.15): Ei(x) = gamma + ln(x) + sum x^k / (k*k!). */
static double ei_series(long double x)
{
    long double sum = EULER + logl(x);
    long double term = 1.0L;

    for (int k = 1; k <= MAXIT; ++k) {
        /* From t_(k-1) to t_k; no factorials or powers are formed. */
        if (k == 1)
            term = x;
        else
            term *= x * (long double)(k - 1) / ((long double)k * k);

        sum += term;
        if (fabsl(term) <= EPS * fabsl(sum))
            break;
    }

    return (double)sum;
}

/* (2.16): e^x/x times an ASYMPTOTIC (not convergent) factorial sum.
   Return 0 if its smallest term is too large for the requested EPS. */
static int ei_asymptotic(long double x, double *result)
{
    long double sum = 1.0L;
    long double term = 1.0L;

    for (int k = 1; k <= MAXIT; ++k) {
        long double next = term * (long double)k / x;

        /* The series ultimately diverges: stop at its smallest term. */
        if (next > term)
            break;

        term = next;
        sum += term;
        if (term == 0.0L)
            break;
    }

    if (term > EPS)
        return 0;

    /* exp(x)/x overflows prematurely near x=710 in double, although
       Ei(x) is still representable. Evaluate the combined exponent. */
    *result = (double)(expl(x - logl(x)) * sum);
    return 1;
}

double Ei(double x)
{
    if (isnan(x) || x <= 0.0)
        return NAN;
    if (isinf(x))
        return INFINITY;

    long double xx = x;
    double approximation;

    /* The textbook suggests |ln(EPS)| as the switching scale.
       Near this boundary an optimal asymptotic truncation may still
       be less accurate than EPS, so fall back to the convergent sum. */
    if (xx >= -logl(EPS) && ei_asymptotic(xx, &approximation))
        return approximation;

    return ei_series(xx);
}

int ei_write_table(FILE *stream, double xmin, double xmax, size_t count)
{
    if (stream == NULL || !isfinite(xmin) || !isfinite(xmax)
        || xmin <= 0.0 || xmax < xmin
        || count == 0 || count > EI_MAX_POINTS)
        return -1;

    for (size_t i = 0; i < count; ++i) {
        double x = xmin;
        if (count > 1) {
            long double fraction = (long double)i / (count - 1);
            long double current = (long double)xmin
                                  + ((long double)xmax - xmin) * fraction;
            x = (i == count - 1) ? xmax : (double)current;
        }

        if (fprintf(stream, "%.17g %.17g\n", x, Ei(x)) < 0)
            return -1;
    }

    return 0;
}
