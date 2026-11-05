#include "linal.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    /* A: 2 x 3, B: 3 x 2, result C: 2 x 2. */
    double **a = alloc2d(2, 3);
    double **b = alloc2d(3, 2);
    double **c = alloc2d(2, 2);
    if (a == NULL || b == NULL || c == NULL) {
        fprintf(stderr, "Error: cannot allocate matrices.\n");
        free2d(a, 2, 3);
        free2d(b, 3, 2);
        free2d(c, 2, 2);
        return EXIT_FAILURE;
    }

    const double left[2][3] = {{1, 2, 3}, {4, 5, 6}};
    const double right[3][2] = {{7, 8}, {9, 10}, {11, 12}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            a[i][j] = left[i][j];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            b[i][j] = right[i][j];

    fill(c, 2, 2, 0.0);
    matrix_prod(a, b, c, 2, 3);
    if (errno != 0) {
        fprintf(stderr, "Error: matrix multiplication failed.\n");
        free2d(a, 2, 3);
        free2d(b, 3, 2);
        free2d(c, 2, 2);
        return EXIT_FAILURE;
    }

    printf("A (2x3):\n");
    show2d(a, 2, 3);
    printf("B (3x2):\n");
    show2d(b, 3, 2);
    printf("C = A*B (2x2):\n");
    show2d(c, 2, 2);

    double u[2] = {3.0, 4.0};
    double v[2] = {1.0, 2.0};
    double length = vabs(u, 2);
    double scalar = scal_prod(u, v, 2);
    if (!isfinite(length) || !isfinite(scalar)) {
        fprintf(stderr, "Error: vector operation failed.\n");
        free2d(a, 2, 3);
        free2d(b, 3, 2);
        free2d(c, 2, 2);
        return EXIT_FAILURE;
    }
    printf("|u| = %.17g\n", length);
    printf("u dot v = %.17g\n", scalar);

    free2d(a, 2, 3);
    free2d(b, 3, 2);
    free2d(c, 2, 2);
    return EXIT_SUCCESS;
}
