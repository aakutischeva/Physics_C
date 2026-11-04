#include "equations.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void print_real_roots(const char *name, double x1, double x2)
{
    if (isnan(x1) || isnan(x2)) {
        printf("%s: real roots are unavailable\n", name);
        return;
    }

    printf("%s: x1 = %.16g, x2 = %.16g\n", name, x1, x2);
}

static void print_complex_root(const char *name, double complex z)
{
    printf("%s = %.16g %+.16gi\n", name, creal(z), cimag(z));
}

int main(void)
{
    int choice;

    printf("Section 2.1: quadratic and cubic equations (C11)\n");
    printf("1 - quadratic: a*x^2 + b*x + c = 0\n");
    printf("2 - cubic: x^3 + a*x^2 + b*x + c = 0\n");
    printf("Choice: ");

    if (scanf("%d", &choice) != 1) {
        fprintf(stderr, "Error: enter 1 or 2.\n");
        return EXIT_FAILURE;
    }

    if (choice == 1) {
        double a, b, c, x1, x2;
        double complex z1, z2;

        printf("Enter a b c (decimal separator is a dot): ");
        if (scanf("%lf %lf %lf", &a, &b, &c) != 3 || a == 0.0
            || !isfinite(a) || !isfinite(b) || !isfinite(c)) {
            fprintf(stderr, "Error: a must be nonzero; all inputs finite.\n");
            return EXIT_FAILURE;
        }

        eq2(a, b, c, &x1, &x2);
        print_real_roots("eq2 (direct formula)", x1, x2);

        eq2_correct(a, b, c, &x1, &x2);
        print_real_roots("eq2_correct (stable)", x1, x2);

        eq2_correct_complex(a, b, c, &z1, &z2);
        print_complex_root("complex x1", z1);
        print_complex_root("complex x2", z2);
    } else if (choice == 2) {
        double a, b, c;
        double complex z1, z2, z3;

        printf("Enter a b c (leading coefficient of x^3 is 1): ");
        if (scanf("%lf %lf %lf", &a, &b, &c) != 3
            || !isfinite(a) || !isfinite(b) || !isfinite(c)) {
            fprintf(stderr, "Error: all inputs must be finite.\n");
            return EXIT_FAILURE;
        }

        eq3(a, b, c, &z1, &z2, &z3);
        print_complex_root("x1", z1);
        print_complex_root("x2", z2);
        print_complex_root("x3", z3);
    } else {
        fprintf(stderr, "Error: unknown menu item.\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
