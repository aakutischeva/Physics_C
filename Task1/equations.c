#include "equations.h"

#include <math.h>
#include <stddef.h>

#ifndef CMPLX
#define CMPLX(real, imag) ((double complex)((real) + (imag) * I))
#endif

void eq2(double a, double b, double c, double *x1, double *x2)
{
    if (x1 == NULL || x2 == NULL)
        return;

    *x1 = NAN;
    *x2 = NAN;

    if (a == 0.0 || !isfinite(a) || !isfinite(b) || !isfinite(c))
        return;

    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0 || isnan(discriminant))
        return;

    double root = sqrt(discriminant);
    *x1 = (-b + root) / (2.0 * a);
    *x2 = (-b - root) / (2.0 * a);
}

void eq2_correct(double a, double b, double c, double *x1, double *x2)
{
    if (x1 == NULL || x2 == NULL)
        return;

    *x1 = NAN;
    *x2 = NAN;

    if (a == 0.0 || !isfinite(a) || !isfinite(b) || !isfinite(c))
        return;

    long double aa = a;
    long double bb = b;
    long double cc = c;
    long double discriminant = bb * bb - 4.0L * aa * cc;

    if (discriminant < 0.0L)
        return;

    if (c == 0.0) {
        *x1 = 0.0;
        *x2 = (double)(-bb / aa);
        return;
    }

    long double root = sqrtl(discriminant);
    long double q = -0.5L * (bb + copysignl(root, bb));

    if (q == 0.0L) {
        /* Defensive handling of a degenerate intermediate result. */
        *x1 = (double)(-bb / (2.0L * aa));
        *x2 = *x1;
        return;
    }

    *x1 = (double)(q / aa);
    *x2 = (double)(cc / q);
}

void eq2_correct_complex(double a, double b, double c,
                         double complex *x1, double complex *x2)
{
    if (x1 == NULL || x2 == NULL)
        return;

    *x1 = CMPLX(NAN, NAN);
    *x2 = CMPLX(NAN, NAN);

    if (a == 0.0 || !isfinite(a) || !isfinite(b) || !isfinite(c))
        return;

    long double aa = a;
    long double bb = b;
    long double cc = c;
    long double discriminant = bb * bb - 4.0L * aa * cc;

    if (discriminant >= 0.0L) {
        double real1, real2;
        eq2_correct(a, b, c, &real1, &real2);
        *x1 = CMPLX(real1, 0.0);
        *x2 = CMPLX(real2, 0.0);
        return;
    }

    long double real = -bb / (2.0L * aa);
    long double imag = sqrtl(-discriminant) / (2.0L * aa);
    *x1 = CMPLX((double)real, (double)imag);
    *x2 = CMPLX((double)real, (double)-imag);
}

void eq3(double a, double b, double c,
         double complex *x1, double complex *x2, double complex *x3)
{
    if (x1 == NULL || x2 == NULL || x3 == NULL)
        return;

    *x1 = CMPLX(NAN, NAN);
    *x2 = CMPLX(NAN, NAN);
    *x3 = CMPLX(NAN, NAN);

    if (!isfinite(a) || !isfinite(b) || !isfinite(c))
        return;

    long double aa = a;
    long double bb = b;
    long double cc = c;
    long double q = (aa * aa - 3.0L * bb) / 9.0L;
    long double r = (2.0L * aa * aa * aa - 9.0L * aa * bb
                     + 27.0L * cc) / 54.0L;
    long double r2 = r * r;
    long double q3 = q * q * q;

    if (r2 < q3) {
        /* Three distinct real roots: formulas (2.7)-(2.10). */
        long double argument = r / sqrtl(q3);
        if (argument > 1.0L)
            argument = 1.0L;
        if (argument < -1.0L)
            argument = -1.0L;

        long double t = acosl(argument);
        long double pi = acosl(-1.0L);
        long double factor = -2.0L * sqrtl(q);
        long double shift = aa / 3.0L;

        *x1 = CMPLX((double)(factor * cosl(t / 3.0L) - shift), 0.0);
        *x2 = CMPLX((double)(factor * cosl((t + 2.0L * pi) / 3.0L)
                               - shift), 0.0);
        *x3 = CMPLX((double)(factor * cosl((t - 2.0L * pi) / 3.0L)
                               - shift), 0.0);
    } else {
        long double delta = r2 - q3;
        long double big_a = -copysignl(cbrtl(fabsl(r) + sqrtl(delta)), r);
        long double big_b = (big_a == 0.0L) ? 0.0L : q / big_a;
        long double real = -0.5L * (big_a + big_b) - aa / 3.0L;
        long double imag = sqrtl(3.0L) * 0.5L * (big_a - big_b);

        *x1 = CMPLX((double)(big_a + big_b - aa / 3.0L), 0.0);
        *x2 = CMPLX((double)real, (double)imag);
        *x3 = CMPLX((double)real, (double)-imag);
    }
}
