#include "equations.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>

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

static void check_triple(double complex a, double complex b, double complex c,
                         double complex x, double complex y, double complex z,
                         double tolerance, const char *message)
{
    double complex actual[3] = {a, b, c};
    double complex expected[3] = {x, y, z};
    int used[3] = {0, 0, 0};
    int all_matched = 1;

    for (int i = 0; i < 3; ++i) {
        int matched = 0;
        for (int j = 0; j < 3; ++j) {
            if (!used[j] && close_complex(actual[j], expected[i], tolerance)) {
                used[j] = 1;
                matched = 1;
                break;
            }
        }
        if (!matched)
            all_matched = 0;
    }

    check(all_matched, message);
}

static void check_cubic_residuals(double a, double b, double c,
                                  double complex x, double complex y,
                                  double complex z)
{
    double complex roots[3] = {x, y, z};
    for (int i = 0; i < 3; ++i) {
        double complex r = roots[i];
        double complex value = ((r + a) * r + b) * r + c;
        double scale = 1.0 + pow(cabs(r), 3.0)
                       + fabs(a) * pow(cabs(r), 2.0)
                       + fabs(b) * cabs(r) + fabs(c);
        check(cabs(value) <= 1e-11 * scale,
              "each cubic root must satisfy the original equation");
    }
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
    RUN_TEST(test_eq3_three_real);
    RUN_TEST(test_eq3_complex_pair);
    RUN_TEST(test_eq3_double_root);
    RUN_TEST(test_eq3_triple_root);
    RUN_TEST(test_eq3_zero_coefficients);
    RUN_TEST(test_eq3_nonfinite_input);

    printf("\nTests: %d; checks: %d; failures: %d\n",
           tests_run, checks, failures);
    return failures == 0 ? 0 : 1;
}
