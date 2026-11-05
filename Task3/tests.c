#include "statistics.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned checks = 0;
static unsigned failures = 0;
static unsigned groups = 0;

#define CHECK(condition, description)                                        \
    do {                                                                     \
        ++checks;                                                            \
        if (!(condition)) {                                                  \
            ++failures;                                                      \
            fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__,     \
                    description);                                            \
        }                                                                    \
    } while (0)

static int close_to(double actual, double expected, double tolerance)
{
    return isfinite(actual) && fabs(actual - expected) <= tolerance;
}

static FILE *sample(const char *text)
{
    FILE *file = tmpfile();
    if (file == NULL)
        return NULL;
    if (fputs(text, file) == EOF) {
        fclose(file);
        return NULL;
    }
    rewind(file);
    return file;
}

static void test_three_values(void)
{
    ++groups;
    FILE *file = sample("1\n2\n3\n");
    CHECK(file != NULL, "temporary file created");
    if (file == NULL)
        return;

    Statistics s;
    CHECK(statistics_from_stream(file, &s) == STATS_OK, "three values accepted");
    CHECK(s.count == 3, "N=3");
    CHECK(close_to(s.mean, 2.0, 1e-14), "mean=2");
    CHECK(close_to(s.variance, 1.0, 1e-14), "sample variance=1");
    CHECK(close_to(s.stddev, 1.0, 1e-14), "standard deviation=1");
    CHECK(close_to(s.skewness, 0.0, 1e-14), "symmetric skewness=0");
    CHECK(close_to(s.excess_kurtosis, -7.0 / 3.0, 1e-14),
          "book's excess kurtosis=-7/3");
    fclose(file);
}

static void test_asymmetry(void)
{
    ++groups;
    FILE *file = sample("0\n0\n0\n10\n");
    CHECK(file != NULL, "right-tail file opened");
    if (file == NULL)
        return;
    Statistics s;
    CHECK(statistics_from_stream(file, &s) == STATS_OK, "right tail accepted");
    CHECK(s.count == 4, "N=4");
    CHECK(close_to(s.mean, 2.5, 1e-13), "mean=2.5");
    CHECK(close_to(s.variance, 25.0, 1e-13), "variance=25");
    CHECK(close_to(s.stddev, 5.0, 1e-13), "sigma=5");
    CHECK(close_to(s.skewness, 0.75, 1e-13), "right tail has positive skew");
    CHECK(close_to(s.excess_kurtosis, -1.6875, 1e-13), "excess=-1.6875");
    fclose(file);

    file = sample("0 10 10 10");
    CHECK(file != NULL, "left-tail file opened");
    if (file == NULL)
        return;
    CHECK(statistics_from_stream(file, &s) == STATS_OK, "left tail accepted");
    CHECK(close_to(s.skewness, -0.75, 1e-13), "left tail has negative skew");
    CHECK(close_to(s.excess_kurtosis, -1.6875, 1e-13),
          "reflection preserves kurtosis");
    fclose(file);
}

static void test_two_values(void)
{
    ++groups;
    FILE *file = sample("1 3");
    CHECK(file != NULL, "two values file opened");
    if (file == NULL)
        return;
    Statistics s;
    CHECK(statistics_from_stream(file, &s) == STATS_OK, "N=2 accepted");
    CHECK(close_to(s.mean, 2.0, 1e-14), "mean=2");
    CHECK(close_to(s.variance, 2.0, 1e-14), "variance uses N-1");
    CHECK(close_to(s.stddev, sqrt(2.0), 1e-14), "standard deviation=sqrt(2)");
    CHECK(close_to(s.skewness, 0.0, 1e-14), "skewness=0");
    CHECK(close_to(s.excess_kurtosis, -2.75, 1e-14), "excess=-2.75");
    fclose(file);
}

static void test_shifted_and_scaled(void)
{
    ++groups;
    FILE *file = sample("1000000000000\n1000000000001\n1000000000002\n");
    CHECK(file != NULL, "shifted input opened");
    if (file == NULL)
        return;
    Statistics s;
    CHECK(statistics_from_stream(file, &s) == STATS_OK, "large offset accepted");
    CHECK(close_to(s.mean, 1000000000001.0, 1e-4), "mean with large offset");
    CHECK(close_to(s.variance, 1.0, 1e-12), "variance with large offset");
    CHECK(close_to(s.skewness, 0.0, 1e-12), "skew with large offset");
    fclose(file);

    file = sample("2\n4\n6\n");
    CHECK(file != NULL, "scaled input opened");
    if (file == NULL)
        return;
    CHECK(statistics_from_stream(file, &s) == STATS_OK, "scaled input accepted");
    CHECK(close_to(s.mean, 4.0, 1e-14), "scaled mean");
    CHECK(close_to(s.variance, 4.0, 1e-14), "variance scales quadratically");
    CHECK(close_to(s.stddev, 2.0, 1e-14), "sigma scales linearly");
    CHECK(close_to(s.excess_kurtosis, -7.0 / 3.0, 1e-14),
          "excess invariant under scaling");
    fclose(file);
}

static void test_degenerate_samples(void)
{
    ++groups;
    Statistics s;
    FILE *file = sample("");
    CHECK(file != NULL, "empty file opened");
    if (file == NULL)
        return;
    CHECK(statistics_from_stream(file, &s) == STATS_EMPTY, "empty file rejected");
    CHECK(s.count == 0, "empty count=0");
    fclose(file);

    file = sample("   \n\t \n");
    CHECK(file != NULL, "whitespace file opened");
    if (file == NULL)
        return;
    CHECK(statistics_from_stream(file, &s) == STATS_EMPTY,
          "whitespace-only file is empty");
    fclose(file);

    file = sample("42\n");
    CHECK(file != NULL, "single value opened");
    if (file == NULL)
        return;
    CHECK(statistics_from_stream(file, &s) == STATS_SINGLE_VALUE,
          "N=1 cannot have sample variance");
    CHECK(close_to(s.mean, 42.0, 0), "mean of single value is defined");
    CHECK(isnan(s.variance) && isnan(s.skewness)
          && isnan(s.excess_kurtosis), "undefined moments are not zero");
    fclose(file);

    file = sample("5 5 5 5");
    CHECK(file != NULL, "constant input opened");
    if (file == NULL)
        return;
    CHECK(statistics_from_stream(file, &s) == STATS_ZERO_VARIANCE,
          "constant data detected");
    CHECK(close_to(s.mean, 5.0, 0), "constant mean=5");
    CHECK(close_to(s.variance, 0.0, 0), "constant variance=0");
    CHECK(close_to(s.stddev, 0.0, 0), "constant sigma=0");
    CHECK(isnan(s.skewness) && isnan(s.excess_kurtosis),
          "standardized moments undefined at sigma=0");
    fclose(file);
}

static void test_bad_input(void)
{
    ++groups;
    Statistics s;
    CHECK(statistics_from_stream(NULL, &s) == STATS_INVALID_ARGUMENT,
          "NULL input rejected");
    FILE *file = sample("1");
    CHECK(file != NULL, "temporary file opened");
    if (file == NULL)
        return;
    CHECK(statistics_from_stream(file, NULL) == STATS_INVALID_ARGUMENT,
          "NULL output rejected");
    fclose(file);

    const char *invalid[] = {"1 2 abc 3", "1 nan 2", "1 inf 2",
                             "1e309", "12junk", "1,5"};
    for (size_t i = 0; i < sizeof invalid / sizeof invalid[0]; ++i) {
        file = sample(invalid[i]);
        CHECK(file != NULL, "invalid input file opened");
        if (file == NULL)
            continue;
        CHECK(statistics_from_stream(file, &s) == STATS_INVALID_DATA,
              "invalid token rejected");
        fclose(file);
    }
}

static void test_numeric_extremes(void)
{
    ++groups;
    FILE *file = sample("1e200\n-1e200\n");
    CHECK(file != NULL, "extreme input file opened");
    if (file == NULL)
        return;
    Statistics s;
    CHECK(statistics_from_stream(file, &s) == STATS_NUMERIC_RANGE,
          "variance larger than DBL_MAX reported rather than infinity");
    fclose(file);

    file = sample("1e-200 2e-200 3e-200");
    CHECK(file != NULL, "tiny variance file opened");
    if (file == NULL)
        return;
    CHECK(statistics_from_stream(file, &s) == STATS_NUMERIC_RANGE,
          "positive variance underflow is detected");
    fclose(file);
}

/* A deterministic normal sample, not dependent on the platform's rand(). */
static uint32_t next_random(uint32_t *state)
{
    *state = *state * UINT32_C(1664525) + UINT32_C(1013904223);
    return *state;
}

static void test_normal_distribution(void)
{
    ++groups;
    FILE *file = tmpfile();
    CHECK(file != NULL, "normal sample file opened");
    if (file == NULL)
        return;

    uint32_t state = UINT32_C(20261105);
    const double two_pi = 6.2831853071795864769;
    for (int i = 0; i < 100000; ++i) {
        double u1 = ((double)next_random(&state) + 0.5) / 4294967296.0;
        double u2 = ((double)next_random(&state) + 0.5) / 4294967296.0;
        double z = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
        if (fprintf(file, "%.17g\n", 5.0 + 2.0 * z) < 0) {
            CHECK(0, "writing normal sample");
            fclose(file);
            return;
        }
    }
    rewind(file);

    Statistics s;
    CHECK(statistics_from_stream(file, &s) == STATS_OK,
          "100000 normally distributed values processed in one pass");
    CHECK(s.count == 100000U, "normal sample size correct");
    CHECK(close_to(s.mean, 5.0, 0.05), "normal mean close to 5");
    CHECK(close_to(s.variance, 4.0, 0.12), "normal sample variance close to 4");
    CHECK(fabs(s.skewness) < 0.09, "normal skew close to zero");
    CHECK(fabs(s.excess_kurtosis) < 0.15, "normal excess close to zero");
    fclose(file);
}

int main(void)
{
    test_three_values();
    test_asymmetry();
    test_two_values();
    test_shifted_and_scaled();
    test_degenerate_samples();
    test_bad_input();
    test_numeric_extremes();
    test_normal_distribution();

    printf("Test groups: %u; checks: %u; failures: %u\n",
           groups, checks, failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
