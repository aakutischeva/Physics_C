#include "equations.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

/* Some older MinGW <complex.h> headers omit the C11 CMPLX macro. */
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

    /* Keep the textbook's naive double-precision calculation unchanged. */
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

    /* Long double reduces intermediate overflow and rounding errors. */
    long double aa = a;
    long double bb = b;
    long double cc = c;

    /* When c=0, q may be zero and c/q would be undefined. */
    if (c == 0.0) {
        *x1 = 0.0;
        *x2 = (double)(-bb / aa);
        return;
    }

    /* Avoid forming a potentially underflowing discriminant when b=0. */
    if (b == 0.0) {
        long double squared_root = -cc / aa;
        if (squared_root >= 0.0L) {
            long double root = sqrtl(squared_root);
            *x1 = (double)root;
            *x2 = (double)-root;
        }
        return;
    }

    long double discriminant = bb * bb - 4.0L * aa * cc;
    if (discriminant < 0.0L)
        return;

    long double root = sqrtl(discriminant);
    long double q = -0.5L * (bb + copysignl(root, bb));

    if (q == 0.0L) {
        /* Do not fabricate a repeated root on arithmetic underflow. */
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

    if (b == 0.0) {
        long double squared_root = -cc / aa;
        if (squared_root >= 0.0L) {
            long double root = sqrtl(squared_root);
            *x1 = CMPLX((double)root, 0.0);
            *x2 = CMPLX((double)-root, 0.0);
        } else {
            long double imag = sqrtl(-squared_root);
            *x1 = CMPLX(0.0, (double)imag);
            *x2 = CMPLX(0.0, (double)-imag);
        }
        return;
    }

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

typedef struct {
    long double re;
    long double im;
} root_ld;

static long double cubic_value(long double x, long double a,
                               long double b, long double c)
{
    return ((x + a) * x + b) * x + c;
}

static long double cubic_derivative(long double x, long double a,
                                    long double b)
{
    return (3.0L * x + 2.0L * a) * x + b;
}

/* Cardano's formulas (2.6)-(2.13) provide the initial roots. */
static void cubic_cardano(long double a, long double b, long double c,
                          root_ld roots[3])
{
    long double q = (a * a - 3.0L * b) / 9.0L;
    long double r = (2.0L * a * a * a - 9.0L * a * b
                     + 27.0L * c) / 54.0L;
    long double r2 = r * r;
    long double q3 = q * q * q;

    if (r2 < q3) {
        long double argument = r / sqrtl(q3);
        if (argument > 1.0L)
            argument = 1.0L;
        if (argument < -1.0L)
            argument = -1.0L;

        long double t = acosl(argument);
        long double pi = acosl(-1.0L);
        long double factor = -2.0L * sqrtl(q);
        long double shift = a / 3.0L;

        roots[0] = (root_ld){factor * cosl(t / 3.0L) - shift, 0.0L};
        roots[1] = (root_ld){factor * cosl((t + 2.0L * pi) / 3.0L)
                             - shift, 0.0L};
        roots[2] = (root_ld){factor * cosl((t - 2.0L * pi) / 3.0L)
                             - shift, 0.0L};
    } else {
        long double delta = r2 - q3;
        long double big_a = -copysignl(cbrtl(fabsl(r) + sqrtl(delta)), r);
        long double big_b = (big_a == 0.0L) ? 0.0L : q / big_a;
        long double real = -0.5L * (big_a + big_b) - a / 3.0L;
        long double imag = sqrtl(3.0L) * 0.5L * (big_a - big_b);

        roots[0] = (root_ld){big_a + big_b - a / 3.0L, 0.0L};
        roots[1] = (root_ld){real, imag};
        roots[2] = (root_ld){real, -imag};
    }
}

/* Newton correction prevents inaccurate deflation of a large root. */
static long double cubic_refine(long double x, long double a,
                                long double b, long double c)
{
    for (int i = 0; i < 40; ++i) {
        long double value = cubic_value(x, a, b, c);
        long double deriv = cubic_derivative(x, a, b);
        if (!isfinite(value) || !isfinite(deriv) || deriv == 0.0L)
            break;

        long double step = value / deriv;
        long double next = x - step;
        if (!isfinite(next) || next == x)
            break;
        x = next;
    }
    return x;
}

/* A cubic has at least one real root; bisection is a safe fallback. */
static long double cubic_bisect(long double a, long double b, long double c)
{
    long double bound = 1.0L + fmaxl(fabsl(a), fmaxl(fabsl(b), fabsl(c)));
    long double left = -bound;
    long double right = bound;

    for (int i = 0; i < 2048; ++i) {
        long double middle = left + (right - left) / 2.0L;
        if (middle == left || middle == right)
            return middle;

        long double value = cubic_value(middle, a, b, c);
        if (value == 0.0L)
            return middle;
        if (value < 0.0L)
            left = middle;
        else
            right = middle;
    }
    return left + (right - left) / 2.0L;
}

/* Divide by (x-real_root) and solve the remaining quadratic stably. */
static void cubic_deflate(long double real_root, long double a,
                          long double b, long double c, root_ld roots[3])
{
    long double linear = a + real_root;
    long double constant = b + real_root * linear;

    /* Vieta's product avoids subtracting two almost equal huge terms. */
    if (real_root != 0.0L) {
        long double from_product = -c / real_root;
        if (isfinite(from_product))
            constant = from_product;

        /* If a and the large root cancel, a+root loses the small
           quadratic coefficient; recover it using b=d-root*linear. */
        if (fabsl(linear) <= 32.0L * LDBL_EPSILON
                              * (fabsl(a) + fabsl(real_root))) {
            long double from_middle = (constant - b) / real_root;
            if (isfinite(from_middle))
                linear = from_middle;
        }
    }
    long double discr = linear * linear - 4.0L * constant;

    /* An exactly repeated root can give a tiny negative discriminant. */
    long double scale = linear * linear + 4.0L * fabsl(constant);
    if (fabsl(discr) <= 16.0L * LDBL_EPSILON * scale)
        discr = 0.0L;

    roots[0] = (root_ld){real_root, 0.0L};
    if (discr < 0.0L) {
        long double real = -linear / 2.0L;
        long double imag = sqrtl(-discr) / 2.0L;
        roots[1] = (root_ld){real, imag};
        roots[2] = (root_ld){real, -imag};
    } else if (constant == 0.0L) {
        roots[1] = (root_ld){0.0L, 0.0L};
        roots[2] = (root_ld){-linear, 0.0L};
    } else {
        long double q = -0.5L * (linear + copysignl(sqrtl(discr), linear));
        roots[1] = (root_ld){q, 0.0L};
        roots[2] = (root_ld){constant / q, 0.0L};
    }
}

/* Relative backward error from substituting each root into the cubic. */
static long double cubic_error(root_ld z, long double a,
                               long double b, long double c)
{
    long double p_re = z.re + a;
    long double p_im = z.im;
    long double t_re = p_re * z.re - p_im * z.im + b;
    long double t_im = p_re * z.im + p_im * z.re;
    p_re = t_re * z.re - t_im * z.im + c;
    p_im = t_re * z.im + t_im * z.re;

    long double magnitude = hypotl(z.re, z.im);
    long double scale = magnitude * magnitude * magnitude
                        + fabsl(a) * magnitude * magnitude
                        + fabsl(b) * magnitude + fabsl(c);
    long double residual = hypotl(p_re, p_im);
    if (scale == 0.0L)
        return residual == 0.0L ? 0.0L : INFINITY;
    return residual / scale;
}

static long double cubic_max_error(const root_ld roots[3], long double a,
                                   long double b, long double c)
{
    long double worst = 0.0L;
    for (int i = 0; i < 3; ++i) {
        long double error = cubic_error(roots[i], a, b, c);
        if (!isfinite(error))
            return INFINITY;
        if (error > worst)
            worst = error;
    }
    return worst;
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

    long double aa = a, bb = b, cc = c;
    root_ld initial[3], result[3], fallback[3];
    cubic_cardano(aa, bb, cc, initial);

    /* For three real roots, prefer the largest-magnitude Cardano root. */
    long double real_root = initial[0].re;
    for (int i = 1; i < 3; ++i) {
        if (initial[i].im == 0.0L && fabsl(initial[i].re) > fabsl(real_root))
            real_root = initial[i].re;
    }

    real_root = cubic_refine(real_root, aa, bb, cc);
    cubic_deflate(real_root, aa, bb, cc, result);

    /* Cardano can lose small roots when coefficients differ greatly. */
    long double score = cubic_max_error(result, aa, bb, cc);
    if (score > 128.0L * LDBL_EPSILON) {
        long double bracketed_root = cubic_bisect(aa, bb, cc);
        bracketed_root = cubic_refine(bracketed_root, aa, bb, cc);
        cubic_deflate(bracketed_root, aa, bb, cc, fallback);
        if (cubic_max_error(fallback, aa, bb, cc) < score) {
            for (int i = 0; i < 3; ++i)
                result[i] = fallback[i];
        }
    }

    *x1 = CMPLX((double)result[0].re, (double)result[0].im);
    *x2 = CMPLX((double)result[1].re, (double)result[1].im);
    *x3 = CMPLX((double)result[2].re, (double)result[2].im);
}
