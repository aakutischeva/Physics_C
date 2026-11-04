#include "equations.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>

static int checks = 0;
static int failures = 0;
static int tests_run = 0;

static void check(int condition, const char *message)
{
    ++checks;
    if (!condition) {
        ++failures;
        fprintf(stderr, "    FAIL: %s\n", message);
    }
}

static int close_real(double actual, double expected,
                      double rel_tol, double abs_tol)
{
    return isfinite(actual)
           && fabs(actual - expected) <= abs_tol + rel_tol * fabs(expected);
}

static int close_complex(double complex actual, double complex expected,
                         double tolerance)
{
    return isfinite(creal(actual)) && isfinite(cimag(actual))
           && cabs(actual - expected) <= tolerance * (1.0 + cabs(expected));
}

static int close_complex_precise(double complex actual, double complex expected,
                                 double rel_tol, double abs_tol)
{
    return isfinite(creal(actual)) && isfinite(cimag(actual))
           && cabs(actual - expected) <= abs_tol + rel_tol * cabs(expected);
}

static void check_pair(double a, double b, double x, double y,
                       double tolerance, const char *message)
{
    int direct = close_real(a, x, tolerance, tolerance)
                 && close_real(b, y, tolerance, tolerance);
    int reverse = close_real(a, y, tolerance, tolerance)
                  && close_real(b, x, tolerance, tolerance);
    check(direct || reverse, message);
}

static void check_complex_pair(double complex a, double complex b,
                               double complex x, double complex y,
                               double tolerance, const char *message)
{
    int direct = close_complex(a, x, tolerance)
                 && close_complex(b, y, tolerance);
    int reverse = close_complex(a, y, tolerance)
                  && close_complex(b, x, tolerance);
    check(direct || reverse, message);
}

static int triple_matches(double complex a, double complex b, double complex c,
                          double complex x, double complex y, double complex z,
                          double rel_tol, double abs_tol)
{
    double complex actual[3] = {a, b, c};
    double complex expected[3] = {x, y, z};
    int used[3] = {0, 0, 0};
    int all_matched = 1;

    for (int i = 0; i < 3; ++i) {
        int matched = 0;
        for (int j = 0; j < 3; ++j) {
            if (!used[j] && close_complex_precise(actual[j], expected[i],
                                                 rel_tol, abs_tol)) {
                used[j] = 1;
                matched = 1;
                break;
            }
        }
        if (!matched)
            all_matched = 0;
    }

    return all_matched;
}

static void check_triple(double complex a, double complex b, double complex c,
                         double complex x, double complex y, double complex z,
                         double tolerance, const char *message)
{
    check(triple_matches(a, b, c, x, y, z, tolerance, tolerance), message);
}

static void check_cubic_residuals(double a, double b, double c,
                                  double complex x, double complex y,
                                  double complex z)
{
    double complex roots[3] = {x, y, z};
    for (int i = 0; i < 3; ++i) {
        double complex r = roots[i];
        long double re = creal(r);
        long double im = cimag(r);
        long double t_re = re + a;
        long double t_im = im;
        long double v_re = t_re * re - t_im * im + b;
        long double v_im = t_re * im + t_im * re;
        t_re = v_re * re - v_im * im + c;
        t_im = v_re * im + v_im * re;
        long double radius = hypotl(re, im);
        long double scale = radius * radius * radius
                            + fabsl(a) * radius * radius
                            + fabsl(b) * radius + fabsl(c);
        long double residual = hypotl(t_re, t_im);
        check(isfinite(t_re) && isfinite(t_im)
              && (scale == 0.0L ? residual == 0.0L
                                : residual <= 1e-11L * scale),
              "each cubic root must satisfy the original equation");
    }
}

static void check_quadratic_residual(double a, double b, double c, double x)
{
    long double xx = x;
    long double value = ((long double)a * xx + b) * xx + c;
    long double scale = fabsl(a) * fabsl(xx) * fabsl(xx)
                        + fabsl(b) * fabsl(xx) + fabsl(c);
    check(isfinite(x) && isfinite(value)
          && (scale == 0.0L ? value == 0.0L : fabsl(value) <= 1e-12L * scale),
          "quadratic root must satisfy original equation");
}

static void test_eq2_two_real(void)
{
    double x, y;
    eq2(1.0, -5.0, 6.0, &x, &y);
    check_pair(x, y, 2.0, 3.0, 1e-13, "eq2: distinct real roots");
}

static void test_eq2_repeated(void)
{
    double x, y;
    eq2(1.0, -4.0, 4.0, &x, &y);
    check_pair(x, y, 2.0, 2.0, 1e-13, "eq2: discriminant zero");
}

static void test_precision_loss(void)
{
    double direct_x, direct_y, stable_x, stable_y;
    eq2(1.0, -1.0, 1e-30, &direct_x, &direct_y);
    eq2_correct(1.0, -1.0, 1e-30, &stable_x, &stable_y);

    check(fmin(fabs(direct_x), fabs(direct_y)) < 1e-20,
          "direct formula loses the tiny root near 1e-30");
    check(close_real(fmin(fabs(stable_x), fabs(stable_y)),
                     1e-30, 1e-12, 1e-40),
          "stable formula recovers the tiny root near 1e-30");
    check(close_real(fmax(stable_x, stable_y), 1.0, 1e-13, 1e-13),
          "stable formula also recovers the large root");
}

static void test_eq2_zero_c(void)
{
    double x, y;
    eq2_correct(1.0, -5.0, 0.0, &x, &y);
    check_pair(x, y, 0.0, 5.0, 1e-13, "eq2_correct: c=0");

    eq2_correct(1.0, 0.0, 0.0, &x, &y);
    check_pair(x, y, 0.0, 0.0, 1e-13,
               "eq2_correct: b=c=0, avoid division 0/0");
}

static void test_eq2_zero_b(void)
{
    double x, y;
    eq2_correct(1.0, 0.0, -4.0, &x, &y);
    check_pair(x, y, -2.0, 2.0, 1e-13, "eq2_correct: b=0");
}

static void test_negative_discriminant(void)
{
    double x, y;
    eq2(1.0, 0.0, 1.0, &x, &y);
    check(isnan(x) && isnan(y), "eq2: no real roots -> NaN");

    eq2_correct(1.0, 0.0, 1.0, &x, &y);
    check(isnan(x) && isnan(y), "eq2_correct: no real roots -> NaN");
}

static void test_complex_quadratic(void)
{
    double complex x, y;
    eq2_correct_complex(1.0, 0.0, 1.0, &x, &y);
    check_complex_pair(x, y, I, -I, 1e-13, "complex quadratic: +/-i");

    eq2_correct_complex(1.0, 2.0, 5.0, &x, &y);
    check_complex_pair(x, y, -1.0 + 2.0 * I, -1.0 - 2.0 * I,
                       1e-13, "complex quadratic: -1 +/- 2i");

    eq2_correct_complex(1.0, -5.0, 6.0, &x, &y);
    check_complex_pair(x, y, 2.0, 3.0, 1e-13,
                       "complex-capable function also handles real roots");
}

static void test_degenerate_quadratic(void)
{
    double x, y;
    double complex z1, z2;
    eq2(0.0, 2.0, 1.0, &x, &y);
    check(isnan(x) && isnan(y), "eq2: a=0 -> invalid quadratic");
    eq2_correct(0.0, 2.0, 1.0, &x, &y);
    check(isnan(x) && isnan(y), "eq2_correct: a=0 -> invalid quadratic");
    eq2_correct_complex(0.0, 2.0, 1.0, &z1, &z2);
    check(isnan(creal(z1)) && isnan(creal(z2)),
          "complex quadratic: a=0 -> invalid quadratic");
}

static void test_large_quadratic_coefficients(void)
{
    double x, y;
    eq2_correct(1e200, -3e200, 2e200, &x, &y);
    check_pair(x, y, 1.0, 2.0, 1e-13,
               "stable formula: large coefficients, finite roots");
    check_quadratic_residual(1e200, -3e200, 2e200, x);
    check_quadratic_residual(1e200, -3e200, 2e200, y);
}

static void test_subnormal_quadratic(void)
{
    double tiny = nextafter(0.0, 1.0);
    double x, y;
    double complex z1, z2;

    eq2_correct(tiny, 0.0, -tiny, &x, &y);
    check_pair(x, y, -1.0, 1.0, 1e-13,
               "stable quadratic: subnormal coefficients, real roots");
    check_quadratic_residual(tiny, 0.0, -tiny, x);
    check_quadratic_residual(tiny, 0.0, -tiny, y);

    eq2_correct_complex(tiny, 0.0, tiny, &z1, &z2);
    check_complex_pair(z1, z2, I, -I, 1e-13,
                       "stable quadratic: subnormal coefficients, complex roots");
}

static void test_nearly_repeated_quadratic(void)
{
    double x, y;
    double c = nextafter(1.0, 0.0);
    double shift = sqrt(1.0 - c);
    eq2_correct(1.0, -2.0, c, &x, &y);
    check_pair(x, y, 1.0 - shift, 1.0 + shift, 1e-12,
               "stable quadratic: small positive discriminant");
    check_quadratic_residual(1.0, -2.0, c, x);
    check_quadratic_residual(1.0, -2.0, c, y);
}

static void test_negative_quadratic_leading_coefficient(void)
{
    double x, y;
    eq2_correct(-1.0, 0.0, 1.0, &x, &y);
    check_pair(x, y, -1.0, 1.0, 1e-13,
               "stable quadratic: negative leading coefficient");
}

static void test_nonfinite_quadratic(void)
{
    double x, y;
    double complex z1, z2;
    eq2(NAN, 2.0, 1.0, &x, &y);
    check(isnan(x) && isnan(y), "eq2: NaN input rejected");
    eq2_correct(INFINITY, 2.0, 1.0, &x, &y);
    check(isnan(x) && isnan(y), "eq2_correct: infinity input rejected");
    eq2_correct_complex(1.0, NAN, 1.0, &z1, &z2);
    check(isnan(creal(z1)) && isnan(creal(z2)),
          "eq2_correct_complex: NaN input rejected");
}

static void test_eq3_three_real(void)
{
    double complex x, y, z;
    eq3(-6.0, 11.0, -6.0, &x, &y, &z);
    check_triple(x, y, z, 1.0, 2.0, 3.0, 1e-12,
                 "eq3: three distinct real roots");
    check_cubic_residuals(-6.0, 11.0, -6.0, x, y, z);
}

static void test_eq3_complex_pair(void)
{
    double complex x, y, z;
    double imag = sqrt(3.0) / 2.0;
    eq3(0.0, 0.0, -1.0, &x, &y, &z);
    check_triple(x, y, z, 1.0, -0.5 + imag * I,
                 -0.5 - imag * I, 1e-12,
                 "eq3: one real root and a conjugate pair");
    check_cubic_residuals(0.0, 0.0, -1.0, x, y, z);
}

static void test_eq3_double_root(void)
{
    double complex x, y, z;
    eq3(0.0, -3.0, 2.0, &x, &y, &z);
    check_triple(x, y, z, -2.0, 1.0, 1.0, 1e-12,
                 "eq3: discriminant boundary, repeated root");
    check_cubic_residuals(0.0, -3.0, 2.0, x, y, z);
}

static void test_eq3_triple_root(void)
{
    double complex x, y, z;
    eq3(-6.0, 12.0, -8.0, &x, &y, &z);
    check_triple(x, y, z, 2.0, 2.0, 2.0, 1e-12,
                 "eq3: triple root, A=0 and B=0");
    check_cubic_residuals(-6.0, 12.0, -8.0, x, y, z);
}

static void test_eq3_zero_coefficients(void)
{
    double complex x, y, z;
    eq3(0.0, 0.0, 0.0, &x, &y, &z);
    check_triple(x, y, z, 0.0, 0.0, 0.0, 1e-12,
                 "eq3: all coefficients zero, root 0 (multiplicity 3)");
    check_cubic_residuals(0.0, 0.0, 0.0, x, y, z);
}

static void test_eq3_nonfinite_input(void)
{
    double complex x, y, z;
    eq3(NAN, 0.0, 0.0, &x, &y, &z);
    check(isnan(creal(x)) && isnan(creal(y)) && isnan(creal(z)),
          "eq3: nonfinite input -> NaN roots");
}

static void test_eq3_double_root_at_zero(void)
{
    double complex x, y, z;
    eq3(1.0, 0.0, 0.0, &x, &y, &z);
    check_triple(x, y, z, -1.0, 0.0, 0.0, 1e-13,
                 "eq3: x^2(x+1), no spurious imaginary part");
    check(cimag(x) == 0.0 && cimag(y) == 0.0 && cimag(z) == 0.0,
          "eq3: exact repeated real roots must remain real");
    check_cubic_residuals(1.0, 0.0, 0.0, x, y, z);
}

static void test_eq3_repeated_negative_root(void)
{
    double complex x, y, z;
    eq3(2.0, 1.0, 0.0, &x, &y, &z);
    check_triple(x, y, z, 0.0, -1.0, -1.0, 1e-13,
                 "eq3: x(x+1)^2, repeated negative root");
    check(cimag(x) == 0.0 && cimag(y) == 0.0 && cimag(z) == 0.0,
          "eq3: repeated negative root must remain real");
    check_cubic_residuals(2.0, 1.0, 0.0, x, y, z);
}

static void test_eq3_separated_roots(void)
{
    double complex x, y, z;
    double large = 1e10;
    eq3(-large, -1.0, large, &x, &y, &z);
    check_triple(x, y, z, -1.0, 1.0, large, 1e-12,
                 "eq3: recover small roots despite one root of size 1e10");
    check_cubic_residuals(-large, -1.0, large, x, y, z);
}

static void test_eq3_extremely_separated_roots(void)
{
    const double large_values[] = {1e20, 1e100, 1e300};
    for (size_t i = 0; i < sizeof(large_values) / sizeof(large_values[0]); ++i) {
        double large = large_values[i];
        double complex x, y, z;
        eq3(-large, -1.0, large, &x, &y, &z);
        check_triple(x, y, z, -1.0, 1.0, large, 1e-12,
                     "eq3: extremely different root magnitudes");
        check_cubic_residuals(-large, -1.0, large, x, y, z);
    }
}

static void test_eq3_two_large_opposite_roots(void)
{
    /* (x - 1e20)(x - 1e100)(x + 1e100) = 0. */
    double complex x, y, z;
    eq3(-1e20, -1e200, 1e220, &x, &y, &z);
    check_triple(x, y, z, 1e20, -1e100, 1e100, 1e-12,
                 "eq3: small real root between two huge opposite roots");
    check_cubic_residuals(-1e20, -1e200, 1e220, x, y, z);
}

static void test_eq3_two_large_imaginary_roots(void)
{
    /* (x - 1e20)(x^2 + 1e200) = 0. */
    double complex x, y, z;
    eq3(-1e20, 1e200, -1e220, &x, &y, &z);
    check_triple(x, y, z, 1e20, 1e100 * I, -1e100 * I, 1e-12,
                 "eq3: small real root with enormous complex roots");
    check_cubic_residuals(-1e20, 1e200, -1e220, x, y, z);
}

static void test_eq3_large_and_two_tiny_real_roots(void)
{
    /* Powers of two make the intended small roots exactly representable.
       (x+2^500)(x+2^-500)^2 is rounded to the supplied coefficients. */
    double huge = ldexp(1.0, 500);
    double tiny = ldexp(1.0, -500);
    double complex x, y, z;

    eq3(huge, 2.0, tiny, &x, &y, &z);
    check(triple_matches(x, y, z, -huge, -tiny, -tiny,
                         1e-12, 0.0),
          "eq3: two tiny real roots despite cancellation in a+root");
    check(cimag(x) == 0.0 && cimag(y) == 0.0 && cimag(z) == 0.0,
          "eq3: tiny real roots must not become imaginary");
    check_cubic_residuals(huge, 2.0, tiny, x, y, z);

    eq3(-huge, 2.0, -tiny, &x, &y, &z);
    check(triple_matches(x, y, z, huge, tiny, tiny,
                         1e-12, 0.0),
          "eq3: large positive root and two tiny positive roots");
    check_cubic_residuals(-huge, 2.0, -tiny, x, y, z);
}

static void test_eq3_large_complex_pair(void)
{
    double complex x, y, z;
    eq3(-1.0, 1e20, -1e20, &x, &y, &z);
    check_triple(x, y, z, 1.0, 1e10 * I, -1e10 * I, 1e-12,
                 "eq3: one small real root and a large imaginary pair");
    check_cubic_residuals(-1.0, 1e20, -1e20, x, y, z);
}

static void test_eq3_tiny_coefficients(void)
{
    double complex x, y, z;
    double mag = cbrt(1e-300);
    double imag = sqrt(3.0) * mag / 2.0;
    eq3(0.0, 0.0, 1e-300, &x, &y, &z);
    check(triple_matches(x, y, z, -mag, mag / 2.0 + imag * I,
                         mag / 2.0 - imag * I, 1e-12, 1e-120),
          "eq3: roots near 1e-100 are accurate relatively, not just absolutely");
    check_cubic_residuals(0.0, 0.0, 1e-300, x, y, z);
}

static void test_eq3_huge_complex_pair(void)
{
    double complex x, y, z;
    eq3(0.0, 1e300, 1.0, &x, &y, &z);
    check(close_real(creal(x), -1e-300, 1e-12, 1e-320)
          && cimag(x) == 0.0,
          "eq3: real root near 1e-300 with huge complex pair");
    check(close_real(fabs(cimag(y)), 1e150, 1e-12, 0.0)
          && close_real(fabs(cimag(z)), 1e150, 1e-12, 0.0),
          "eq3: complex roots with imaginary magnitude 1e150");
    check_cubic_residuals(0.0, 1e300, 1.0, x, y, z);
}

static uint32_t rng_state = UINT32_C(0x12345678);

static double sample_coefficient(void)
{
    rng_state = UINT32_C(1664525) * rng_state + UINT32_C(1013904223);
    return 10.0 * ((double)rng_state / (double)UINT32_MAX) - 5.0;
}

static void test_eq3_deterministic_random_cases(void)
{
    int all_real_ok = 1;
    int all_complex_ok = 1;
    for (int i = 0; i < 500; ++i) {
        double u = sample_coefficient();
        double v = sample_coefficient();
        double w = sample_coefficient();
        if (fabs(u - v) >= 0.2 && fabs(u - w) >= 0.2
            && fabs(v - w) >= 0.2) {
            double a = -(u + v + w);
            double b = u * v + u * w + v * w;
            double c = -u * v * w;
            double complex x, y, z;
            eq3(a, b, c, &x, &y, &z);
            if (!triple_matches(x, y, z, u, v, w, 1e-8, 1e-8))
                all_real_ok = 0;
        }

        double real = sample_coefficient();
        double part = sample_coefficient();
        double imag = 0.5 + fabs(sample_coefficient());
        double a = -(real + 2.0 * part);
        double b = 2.0 * real * part + part * part + imag * imag;
        double c = -real * (part * part + imag * imag);
        double complex x, y, z;
        eq3(a, b, c, &x, &y, &z);
        if (!triple_matches(x, y, z, real, part + imag * I,
                            part - imag * I, 1e-8, 1e-8))
            all_complex_ok = 0;
    }
    check(all_real_ok, "eq3: 500 reproducible three-real-root cases");
    check(all_complex_ok, "eq3: 500 reproducible conjugate-pair cases");
}

#define RUN_TEST(function) do {                          \
    int before = failures;                               \
    function();                                          \
    ++tests_run;                                         \
    printf("[%s] %s\n", failures == before ? "PASS" : "FAIL", \
           #function);                                  \
} while (0)

int main(void)
{
    RUN_TEST(test_eq2_two_real);
    RUN_TEST(test_eq2_repeated);
    RUN_TEST(test_precision_loss);
    RUN_TEST(test_eq2_zero_c);
    RUN_TEST(test_eq2_zero_b);
    RUN_TEST(test_negative_discriminant);
    RUN_TEST(test_complex_quadratic);
    RUN_TEST(test_degenerate_quadratic);
    RUN_TEST(test_large_quadratic_coefficients);
    RUN_TEST(test_subnormal_quadratic);
    RUN_TEST(test_nearly_repeated_quadratic);
    RUN_TEST(test_negative_quadratic_leading_coefficient);
    RUN_TEST(test_nonfinite_quadratic);
    RUN_TEST(test_eq3_three_real);
    RUN_TEST(test_eq3_complex_pair);
    RUN_TEST(test_eq3_double_root);
    RUN_TEST(test_eq3_triple_root);
    RUN_TEST(test_eq3_zero_coefficients);
    RUN_TEST(test_eq3_nonfinite_input);
    RUN_TEST(test_eq3_double_root_at_zero);
    RUN_TEST(test_eq3_repeated_negative_root);
    RUN_TEST(test_eq3_separated_roots);
    RUN_TEST(test_eq3_extremely_separated_roots);
    RUN_TEST(test_eq3_two_large_opposite_roots);
    RUN_TEST(test_eq3_two_large_imaginary_roots);
    RUN_TEST(test_eq3_large_and_two_tiny_real_roots);
    RUN_TEST(test_eq3_large_complex_pair);
    RUN_TEST(test_eq3_tiny_coefficients);
    RUN_TEST(test_eq3_huge_complex_pair);
    RUN_TEST(test_eq3_deterministic_random_cases);

    printf("\nTests: %d; checks: %d; failures: %d\n",
           tests_run, checks, failures);
    return failures == 0 ? 0 : 1;
}
