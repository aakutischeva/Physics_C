#include "integrals.h"

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned checks = 0;
static unsigned failures = 0;
static unsigned groups = 0;

#define CHECK(test, description)                                           \
    do {                                                                   \
        ++checks;                                                          \
        if (!(test)) {                                                     \
            ++failures;                                                    \
            fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__,    \
                    description);                                          \
        }                                                                  \
    } while (0)

static int near(double actual, double expected, double tol)
{
    return isfinite(actual) && fabs(actual - expected) <= tol;
}

static double one(double x) { (void)x; return 1.0; }
static double zero(double x) { (void)x; return 0.0; }
static double line(double x) { return x; }
static double square(double x) { return x * x; }
static double cubic(double x) { return x * x * x; }
static double fourth(double x) { double t = x * x; return t * t; }
static double inverse_quadratic(double x) { return 1.0 / (1.0 + x * x); }
static double pole(double x) { return 1.0 / x; }
static double not_finite(double x) { (void)x; return NAN; }
static double shifted(double x) { return x - 1e12; }
static double tiny(double x) { (void)x; return nextafter(0.0, 1.0); }

static unsigned int evaluation_count = 0U;
static double counted_one(double x)
{
    (void)x;
    ++evaluation_count;
    return 1.0;
}

static void test_trapezoids(void)
{
    ++groups;
    CHECK(near(trapz(one, 0, 1, 1), 1, 1e-15), "constant exact, n=1");
    CHECK(near(trapz(line, 0, 1, 1), 0.5, 1e-15), "linear exact, n=1");
    CHECK(near(trapz(line, -5, 3, 7), -8, 1e-14), "linear with mixed signs");
    CHECK(near(trapz(square, 0, 1, 2), 0.375, 1e-15), "quadratic trapezoids n=2");
    CHECK(near(trapz(line, 1, 0, 10), -0.5, 1e-14), "reversed bounds");
    CHECK(near(trapz(zero, 0, 1, 10), 0.0, 0.0), "zero function");
    CHECK(trapz(square, 2, 2, 1) == 0.0, "zero-width interval");
    CHECK(near(trapz(shifted, 1e12, 1e12 + 1, 100), 0.5, 1e-4),
          "large common offset");
}

static void test_formula_223(void)
{
    ++groups;
    CHECK(near(simpson(one, 0, 1, 5), 1, 1e-14), "constant exact, n=5");
    CHECK(near(simpson(line, 0, 1, 5), 0.5, 1e-14), "linear exact, n=5");
    CHECK(near(simpson(square, 0, 1, 5), 1.0/3.0, 1e-14),
          "quadratic exact, n=5");
    CHECK(near(simpson(cubic, 0, 1, 5), 0.25, 1e-14),
          "cubic exact, n=5");
    CHECK(near(simpson(cubic, -2, 3, 7), 65.0/4.0, 1e-12),
          "cubic exact with arbitrary bounds");
    CHECK(near(simpson(square, 1, 0, 5), -1.0/3.0, 1e-14),
          "reversed bounds");
    CHECK(simpson(one, 9, 9, 5) == 0.0, "zero-width interval");
    double err10 = fabs(simpson(fourth, 0, 1, 10) - 0.2);
    double err20 = fabs(simpson(fourth, 0, 1, 20) - 0.2);
    CHECK(err10 > 0.0 && err20 < err10 / 10.0,
          "fourth-order convergence on quartic polynomial");
}

static void test_adaptive(void)
{
    ++groups;
    const double pi = acos(-1.0);
    CHECK(near(trapz_adapt(square, 0, 1, 1e-9), 1.0/3.0, 1e-8),
          "adaptive quadratic");
    CHECK(near(trapz_adapt(sin, 0, pi, 1e-9), 2.0, 1e-8),
          "adaptive sine");
    CHECK(near(trapz_adapt(inverse_quadratic, 0, 1, 1e-9), pi/4, 1e-8),
          "adaptive arctangent integral");
    CHECK(near(trapz_adapt(square, 1, 0, 1e-9), -1.0/3.0, 1e-8),
          "adaptive reversed bounds, negative result");
    CHECK(near(trapz_adapt(line, -1, 1, 1e-9), 0.0, 1e-14),
          "adaptive zero integral");
    CHECK(trapz_adapt(square, 7, 7, 1e-8) == 0.0,
          "adaptive zero-width interval");

    evaluation_count = 0U;
    CHECK(near(trapz_adapt(counted_one, 0, 1, 1e-9), 1.0, 1e-14),
          "adaptive constant value");
    CHECK(evaluation_count == 3U,
          "adaptive reuses two endpoints and adds only one midpoint");
}

static void test_adaptive_small_integrals(void)
{
    ++groups;
    /* For x^2 the answer scales as c^3/3.  The stopping criterion
       must preserve relative accuracy even when the integral is tiny. */
    const double bounds[] = {1.0, 0.01, 1e-6, 1e-9, 1e-12};
    for (size_t i = 0; i < sizeof bounds / sizeof bounds[0]; ++i) {
        double c = bounds[i];
        double exact = (c * c * c) / 3.0;
        double actual = trapz_adapt(square, 0.0, c, 1e-6);
        CHECK(isfinite(actual) && fabs(actual / exact - 1.0) < 1e-6,
              "small positive integrals retain relative accuracy");

        double reversed = trapz_adapt(square, c, 0.0, 1e-6);
        CHECK(isfinite(reversed) && fabs(reversed / -exact - 1.0) < 1e-6,
              "reversed small integrals retain relative accuracy");
    }

    CHECK(near(trapz_adapt(line, -1.0, 1.0, 1e-6), 0.0, 1e-14),
          "zero integral is accepted despite a nonzero absolute integral");
    CHECK(trapz_adapt(zero, 0.0, 1.0, 1e-6) == 0.0,
          "identically zero function converges with zero tolerance scale");
}

static void test_invalid_arguments(void)
{
    ++groups;
    CHECK(isnan(trapz(NULL, 0, 1, 5)) && errno == EINVAL,
          "NULL function pointer rejected");
    CHECK(isnan(trapz(line, 0, 1, 0)) && errno == EINVAL,
          "n=0 rejected");
    CHECK(isnan(trapz(line, 0, 1, -1)) && errno == EINVAL,
          "negative n rejected");
    CHECK(isnan(simpson(line, 0, 1, 4)) && errno == EINVAL,
          "formula (2.23) requires n>=5");
    CHECK(isnan(simpson(line, 0, 1, -100)) && errno == EINVAL,
          "negative n for formula (2.23)");
    CHECK(isnan(trapz(line, NAN, 1, 5)) && errno == EINVAL,
          "NaN integration bound rejected");
    CHECK(isnan(trapz(line, 0, INFINITY, 5)) && errno == EINVAL,
          "infinite integration bound rejected");
    CHECK(isnan(trapz_adapt(line, 0, 1, 0)) && errno == EINVAL,
          "zero eps rejected");
    CHECK(isnan(trapz_adapt(line, 0, 1, -0.1)) && errno == EINVAL,
          "negative eps rejected");
    CHECK(isnan(trapz_adapt(line, 0, 1, 1)) && errno == EINVAL,
          "eps>=1 rejected");
    CHECK(isnan(trapz_adapt(line, 0, 1, NAN)) && errno == EINVAL,
          "NaN eps rejected");
    CHECK(isnan(trapz_adapt(NULL, 0, 1, 1e-3)) && errno == EINVAL,
          "NULL adaptive function rejected");
}

static void test_nonfinite_and_range(void)
{
    ++groups;
    CHECK(isnan(trapz(pole, 0, 1, 10)) && errno == EDOM,
          "singular endpoint rejected");
    CHECK(isnan(simpson(not_finite, 0, 1, 5)) && errno == EDOM,
          "non-finite integrand rejected");
    CHECK(isnan(trapz_adapt(pole, 0, 1, 1e-6)) && errno == EDOM,
          "singular adaptive endpoint rejected");
    CHECK(isnan(trapz(one, -DBL_MAX, DBL_MAX, 10)) && errno == ERANGE,
          "unrepresentable integral rejected");
    double smallest = nextafter(0.0, 1.0);
    CHECK(trapz(one, 0.0, smallest, 1) == smallest,
          "smallest subnormal integral is preserved");
    CHECK(isnan(trapz(tiny, 0.0, 0.5, 1)) && errno == ERANGE,
          "nonzero integral underflow to zero is rejected");
}

int main(void)
{
    test_trapezoids();
    test_formula_223();
    test_adaptive();
    test_adaptive_small_integrals();
    test_invalid_arguments();
    test_nonfinite_and_range();

    printf("Test groups: %u; checks: %u; failures: %u\n",
           groups, checks, failures);
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
