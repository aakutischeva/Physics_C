#include "integrals.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double constant(double x) { (void)x; return 1.0; }
static double linear(double x) { return x; }
static double square(double x) { return x * x; }
static double cubic(double x) { return x * x * x; }
static double inverse_quadratic(double x) { return 1.0 / (1.0 + x * x); }

typedef struct {
    const char *name;
    finteg function;
} ExampleFunction;

static const ExampleFunction functions[] = {
    {"one", constant}, {"x", linear}, {"x2", square},
    {"x3", cubic}, {"sin", sin}, {"cos", cos},
    {"exp", exp}, {"invquad", inverse_quadratic}
};

static int parse_double(const char *text, double *value)
{
    char *end;
    errno = 0;
    double parsed = strtod(text, &end);
    if (text == end || *end != '\0' || !isfinite(parsed)
        || (errno == ERANGE && parsed == 0.0))
        return 0;
    *value = parsed;
    return 1;
}

static int parse_count(const char *text, int *value)
{
    char *end;
    errno = 0;
    long parsed = strtol(text, &end, 10);
    if (text == end || *end != '\0' || errno == ERANGE
        || parsed < 5L || parsed > INT_MAX)
        return 0;
    *value = (int)parsed;
    return 1;
}

static void usage(const char *program)
{
    fprintf(stderr, "Usage: %s function a b n eps\n", program);
    fprintf(stderr, "Functions: one x x2 x3 sin cos exp invquad\n");
    fprintf(stderr, "Requirements: n >= 5, 0 < eps < 1.\n");
    fprintf(stderr, "Example: %s x2 0 1 100 1e-9\n", program);
}

int main(int argc, char **argv)
{
    if (argc != 6) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    finteg f = NULL;
    for (size_t i = 0; i < sizeof functions / sizeof functions[0]; ++i) {
        if (strcmp(argv[1], functions[i].name) == 0) {
            f = functions[i].function;
            break;
        }
    }

    double a, b, eps;
    int n;
    if (f == NULL || !parse_double(argv[2], &a)
        || !parse_double(argv[3], &b) || !parse_count(argv[4], &n)
        || !parse_double(argv[5], &eps) || eps <= 0.0 || eps >= 1.0) {
        fprintf(stderr, "Error: invalid command-line argument.\n");
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    double trapezoids = trapz(f, a, b, n);
    if (!isfinite(trapezoids)) {
        fprintf(stderr, "Error: trapezoidal integration failed.\n");
        return EXIT_FAILURE;
    }

    double formula_223 = simpson(f, a, b, n);
    if (!isfinite(formula_223)) {
        fprintf(stderr, "Error: formula (2.23) failed.\n");
        return EXIT_FAILURE;
    }

    double adaptive = trapz_adapt(f, a, b, eps);
    if (!isfinite(adaptive)) {
        fprintf(stderr, "Error: adaptive integration failed or did not converge.\n");
        return EXIT_FAILURE;
    }

    printf("Function: %s\n", argv[1]);
    printf("Interval: [%.17g, %.17g]\n", a, b);
    printf("Intervals: %d; tolerance: %.6g\n", n, eps);
    printf("Trapezoidal (2.21): %.17g\n", trapezoids);
    printf("Formula (2.23):    %.17g\n", formula_223);
    printf("Adaptive:          %.17g\n", adaptive);

    return EXIT_SUCCESS;
}
