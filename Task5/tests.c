#include "linal.h"

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned groups = 0U;
static unsigned checks = 0U;
static unsigned failures = 0U;

#define CHECK(test, message)                                                \
    do {                                                                    \
        ++checks;                                                           \
        if (!(test)) {                                                      \
            ++failures;                                                     \
            fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__,     \
                    message);                                               \
        }                                                                   \
    } while (0)

static int close_to(double x, double expected, double tol)
{
    return isfinite(x) && fabs(x - expected) <= tol;
}

static void test_allocation_and_fill(void)
{
    ++groups;
    double **a = alloc2d(2, 3);
    CHECK(a != NULL, "2x3 allocation");
    if (a != NULL) {
        int zero = 1;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 3; ++j)
                zero = zero && a[i][j] == 0.0;
        CHECK(zero, "matrix starts zero-initialized");

        fill(a, 2, 3, 2.5);
        int filled = 1;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 3; ++j)
                filled = filled && a[i][j] == 2.5;
        CHECK(filled, "fill changes every element");

        fill(a, 2, 3, NAN);
        CHECK(errno == EDOM, "fill rejects NaN");
        CHECK(a[0][0] == 2.5 && a[1][2] == 2.5,
              "invalid fill does not modify matrix");
        free2d(a, 2, 3);
    }

    CHECK(alloc2d(0, 3) == NULL && errno == EINVAL,
          "zero row count rejected");
    CHECK(alloc2d(2, -1) == NULL && errno == EINVAL,
          "negative column count rejected");
    free2d(NULL, 0, 0);
    CHECK(1, "free2d(NULL) is a no-op");
}

static void test_vector_norm(void)
{
    ++groups;
    double a[2] = {3.0, 4.0};
    CHECK(close_to(vabs(a, 2), 5.0, 1e-15), "3-4-5 vector");

    double small[2] = {1e-300, 1e-300};
    CHECK(close_to(vabs(small, 2) / 1e-300, sqrt(2.0), 1e-14),
          "tiny components do not vanish when squared");

    double large[2] = {1e308, 1e308};
    CHECK(close_to(vabs(large, 2) / 1e308, sqrt(2.0), 1e-14),
          "large components do not overflow when squared");

    double overflow[2] = {DBL_MAX, DBL_MAX};
    CHECK(isnan(vabs(overflow, 2)) && errno == ERANGE,
          "out-of-range norm reported");

    double nan_vec[2] = {1.0, NAN};
    CHECK(isnan(vabs(nan_vec, 2)) && errno == EDOM,
          "NaN component rejected");
    CHECK(isnan(vabs(NULL, 2)) && errno == EINVAL,
          "NULL vector rejected");
    CHECK(isnan(vabs(a, 0)) && errno == EINVAL,
          "zero vector length rejected");

    double zeros[3] = {0, 0, 0};
    CHECK(vabs(zeros, 3) == 0.0, "zero vector norm");
}

static void test_scalar_product(void)
{
    ++groups;
    double a[3] = {1, 2, 3};
    double b[3] = {4, 5, 6};
    CHECK(scal_prod(a, b, 3) == 32.0, "three-component dot product");
    CHECK(scal_prod(a, a, 3) == 14.0, "vector squared norm");

    double cancel_a[3] = {1e16, 1, -1e16};
    double ones[3] = {1, 1, 1};
    CHECK(scal_prod(cancel_a, ones, 3) == 1.0,
          "compensated cancellation preserves small term");

    double huge[1] = {1e200};
    CHECK(isnan(scal_prod(huge, huge, 1)) && errno == ERANGE,
          "dot product overflow rejected");

    double nan_a[2] = {1.0, NAN};
    CHECK(isnan(scal_prod(nan_a, b, 2)) && errno == EDOM,
          "non-finite dot input rejected");
    CHECK(isnan(scal_prod(a, NULL, 3)) && errno == EINVAL,
          "NULL second vector rejected");
    CHECK(isnan(scal_prod(a, b, -1)) && errno == EINVAL,
          "negative length rejected");
}

static void test_matrix_product(void)
{
    ++groups;
    double **a = alloc2d(2, 3);
    double **b = alloc2d(3, 2);
    double **c = alloc2d(2, 2);
    CHECK(a != NULL && b != NULL && c != NULL,
          "rectangular product matrices allocated");
    if (a == NULL || b == NULL || c == NULL) {
        free2d(a, 2, 3); free2d(b, 3, 2); free2d(c, 2, 2);
        return;
    }

    const double av[2][3] = {{1, 2, 3}, {4, 5, 6}};
    const double bv[3][2] = {{7, 8}, {9, 10}, {11, 12}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            a[i][j] = av[i][j];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            b[i][j] = bv[i][j];

    matrix_prod(a, b, c, 2, 3);
    CHECK(errno == 0, "matrix_prod succeeds");
    CHECK(c[0][0] == 58 && c[0][1] == 64,
          "first result row correct");
    CHECK(c[1][0] == 139 && c[1][1] == 154,
          "second result row correct");

    fill(c, 2, 2, 42.0);
    matrix_prod(NULL, b, c, 2, 3);
    CHECK(errno == EINVAL, "NULL matrix rejected");
    CHECK(c[0][0] == 42.0 && c[1][1] == 42.0,
          "output unchanged after invalid input");

    a[1][1] = NAN;
    matrix_prod(a, b, c, 2, 3);
    CHECK(errno == EDOM, "non-finite matrix element rejected");
    CHECK(c[0][0] == 42.0 && c[1][1] == 42.0,
          "numeric error does not partially overwrite output");

    free2d(a, 2, 3);
    free2d(b, 3, 2);
    free2d(c, 2, 2);
}

static void test_general_rectangular_product(void)
{
    ++groups;
    double **a = alloc2d(2, 3);
    double **b = alloc2d(3, 4);
    double **c = alloc2d(2, 4);
    CHECK(a != NULL && b != NULL && c != NULL,
          "general rectangular matrices allocated");
    if (a == NULL || b == NULL || c == NULL) {
        free2d(a, 2, 3); free2d(b, 3, 4); free2d(c, 2, 4);
        return;
    }

    const double av[2][3] = {{1, 2, 3}, {4, 5, 6}};
    const double bv[3][4] = {{1, 0, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 1, 1}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            a[i][j] = av[i][j];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            b[i][j] = bv[i][j];

    matrix_prod_rect(a, b, c, 2, 3, 4);
    CHECK(errno == 0, "general rectangular multiplication succeeds");
    int correct = 1;
    const double expected[2][4] = {{1, 2, 3, 3}, {4, 5, 6, 6}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 4; ++j)
            correct = correct && c[i][j] == expected[i][j];
    CHECK(correct, "2x3 times 3x4 gives 2x4");

    free2d(a, 2, 3);
    free2d(b, 3, 4);
    free2d(c, 2, 4);
}

static void test_in_place_product(void)
{
    ++groups;
    double **a = alloc2d(2, 2);
    double **b = alloc2d(2, 2);
    CHECK(a != NULL && b != NULL, "square matrices allocated");
    if (a == NULL || b == NULL) {
        free2d(a, 2, 2); free2d(b, 2, 2);
        return;
    }

    a[0][0]=1; a[0][1]=2; a[1][0]=3; a[1][1]=4;
    b[0][0]=5; b[0][1]=6; b[1][0]=7; b[1][1]=8;
    matrix_prod(a, b, a, 2, 2);
    CHECK(errno == 0, "in-place multiplication succeeds");
    CHECK(a[0][0]==19 && a[0][1]==22 && a[1][0]==43 && a[1][1]==50,
          "in-place output is correct");

    free2d(a, 2, 2);
    free2d(b, 2, 2);
}

static void test_matrix_edge_cases(void)
{
    ++groups;
    double **a = alloc2d(1, 1);
    double **b = alloc2d(1, 1);
    double **c = alloc2d(1, 1);
    CHECK(a != NULL && b != NULL && c != NULL,
          "1x1 matrices allocated");
    if (a == NULL || b == NULL || c == NULL) {
        free2d(a, 1, 1); free2d(b, 1, 1); free2d(c, 1, 1);
        return;
    }

    a[0][0] = -3.0;
    b[0][0] = 2.0;
    matrix_prod(a, b, c, 1, 1);
    CHECK(errno == 0 && c[0][0] == -6.0, "1x1 multiplication");

    b[0][0] = 0.0;
    matrix_prod(a, b, c, 1, 1);
    CHECK(errno == 0 && c[0][0] == 0.0, "multiplication by zero");

    c[0][0] = 17.0;
    matrix_prod_rect(a, b, c, 1, 0, 1);
    CHECK(errno == EINVAL && c[0][0] == 17.0,
          "invalid inner dimension preserves output");

    a[0][0] = 1e200;
    b[0][0] = 1e200;
    matrix_prod(a, b, c, 1, 1);
    CHECK(errno == ERANGE && c[0][0] == 17.0,
          "overflow preserves output");

    a[0][0] = 1e-200;
    b[0][0] = 1e-200;
    matrix_prod(a, b, c, 1, 1);
    CHECK(errno == ERANGE && c[0][0] == 17.0,
          "nonzero underflow preserves output");

    b[0][0] = INFINITY;
    matrix_prod(a, b, c, 1, 1);
    CHECK(errno == EDOM && c[0][0] == 17.0,
          "infinite element rejected without output change");

    show2d(NULL, 1, 1);
    CHECK(errno == EINVAL, "show2d rejects NULL");
    show2d(c, 1, 0);
    CHECK(errno == EINVAL, "show2d rejects invalid dimensions");

    free2d(a, 1, 1);
    free2d(b, 1, 1);
    free2d(c, 1, 1);
}

int main(void)
{
    test_allocation_and_fill();
    test_vector_norm();
    test_scalar_product();
    test_matrix_product();
    test_general_rectangular_product();
    test_in_place_product();
    test_matrix_edge_cases();

    printf("Test groups: %u; checks: %u; failures: %u\n",
           groups, checks, failures);
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
