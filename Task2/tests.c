#include "ei.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned checks = 0;
static unsigned failures = 0;
static unsigned cases = 0;

#define CHECK(condition, description)                                          \
    do {                                                                       \
        ++checks;                                                              \
        if (!(condition)) {                                                    \
            ++failures;                                                        \
            fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__,       \
                    description);                                              \
        }                                                                      \
    } while (0)

static int near(double actual, double expected, double relative, double absolute)
{
    return isfinite(actual)
        && fabs(actual - expected) <= absolute + relative * fabs(expected);
}

static void test_reference_values(void)
{
    /* Independent reference values: scipy.special.expi (double precision). */
    static const struct {
        double x;
        double expected;
    } examples[] = {
        {1e-12, -27.053805451026015},
        {1e-6, -13.23829389306249},
        {0.1, -1.6228128139692766},
        {0.3725, -2.887418318875401e-5},
        {1.0, 1.895117816355937},
        {5.0, 40.185275355803171},
        {10.0, 2492.2289762418773},
        {16.0, 595560.99867083703},
        {16.11809565095832, 664918.40504856897},
        {18.0, 3877904.3305974416},
        {19.0, 9950907.2510468438},
        {20.0, 25615652.664056595},
        {22.0, 171144671.30036366},
        {30.0, 368973209407.2738},
        {40.0, 6039718263611238.0},
        {100.0, 2.7155527448538798e41},
        {700.0, 1.4509787360525605e301}
    };

    ++cases;
    for (size_t i = 0; i < sizeof(examples) / sizeof(examples[0]); ++i) {
        double value = Ei(examples[i].x);
        if (!near(value, examples[i].expected, 1e-7, 1e-12)) {
            fprintf(stderr, "  x=%.17g: got %.17g; expected %.17g\n",
                    examples[i].x, value, examples[i].expected);
        }
        CHECK(near(value, examples[i].expected, 1e-7, 1e-12),
              "Ei(x) agrees with an independent reference");
    }
}

static void test_zero_and_domain(void)
{
    ++cases;
    CHECK(isnan(Ei(0.0)), "Ei(0) is outside the domain");
    CHECK(isnan(Ei(-1.0)), "negative x is outside the assignment domain");
    CHECK(isnan(Ei(NAN)), "NaN input returns NaN");
    CHECK(isinf(Ei(INFINITY)) && Ei(INFINITY) > 0,
          "Ei(+infinity) returns +infinity");
    CHECK(isfinite(Ei(DBL_TRUE_MIN)), "smallest positive double is accepted");
    CHECK(near(Ei(DBL_TRUE_MIN), 0.57721566490153286061 + log(DBL_TRUE_MIN),
               1e-12, 1e-10), "tiny x follows gamma + log(x)");
    CHECK(isfinite(Ei(710.0)) && Ei(710.0) > 1e305,
          "Ei(710) does not overflow prematurely");
    CHECK(isinf(Ei(717.0)) && Ei(717.0) > 0,
          "Ei(717) exceeds double range");
}

static void test_close_to_zero_of_ei(void)
{
    ++cases;
    CHECK(fabs(Ei(0.3725074107813666)) <= 1e-12,
          "near a zero of Ei, use absolute rather than relative tolerance");
}

static void test_table_three_points(void)
{
    ++cases;
    FILE *file = tmpfile();
    CHECK(file != NULL, "temporary stream opened");
    if (file == NULL)
        return;

    CHECK(ei_write_table(file, 1.0, 3.0, 3) == 0, "three-point table written");
    rewind(file);

    char line[256];
    for (int i = 0; i < 3; ++i) {
        double x, value;
        char extra;
        CHECK(fgets(line, sizeof line, file) != NULL, "data row exists");
        if (feof(file))
            break;
        CHECK(sscanf(line, "%lf %lf %c", &x, &value, &extra) == 2,
              "each row has exactly two numeric columns");
        CHECK(near(x, (double)(i + 1), 0, 1e-15),
              "x grid includes endpoints and has equal spacing");
        CHECK(near(value, Ei((double)(i + 1)), 1e-14, 1e-14),
              "table value equals Ei(x)");
    }
    CHECK(fgets(line, sizeof line, file) == NULL, "no extra data rows");
    CHECK(fclose(file) == 0, "temporary stream closed");
}

static void test_single_point_and_equal_endpoints(void)
{
    ++cases;
    FILE *file = tmpfile();
    CHECK(file != NULL, "temporary stream opened for single point");
    if (file == NULL)
        return;

    CHECK(ei_write_table(file, 2.0, 5.0, 1) == 0,
          "count=1 writes xmin without division by zero");
    rewind(file);
    double x, value;
    CHECK(fscanf(file, "%lf %lf", &x, &value) == 2, "single row parsed");
    CHECK(near(x, 2.0, 0, 1e-15), "single row uses xmin");
    CHECK(near(value, Ei(2.0), 1e-14, 1e-14), "single value correct");
    CHECK(fscanf(file, "%lf %lf", &x, &value) == EOF, "exactly one row");
    fclose(file);

    file = tmpfile();
    CHECK(file != NULL, "temporary stream opened for equal endpoints");
    if (file == NULL)
        return;
    CHECK(ei_write_table(file, 2.0, 2.0, 3) == 0,
          "equal endpoints allowed");
    rewind(file);
    for (int i = 0; i < 3; ++i) {
        CHECK(fscanf(file, "%lf %lf", &x, &value) == 2, "repeated row parsed");
        CHECK(near(x, 2.0, 0, 1e-15), "repeated x is correct");
    }
    CHECK(fscanf(file, "%lf %lf", &x, &value) == EOF,
          "equal endpoints produce requested row count");
    fclose(file);
}

static void test_invalid_tables(void)
{
    ++cases;
    FILE *file = tmpfile();
    CHECK(file != NULL, "temporary stream opened for invalid arguments");
    if (file == NULL)
        return;

    CHECK(ei_write_table(NULL, 1.0, 2.0, 5) == -1, "NULL stream rejected");
    CHECK(ei_write_table(file, 0.0, 2.0, 5) == -1, "xmin=0 rejected");
    CHECK(ei_write_table(file, -1.0, 2.0, 5) == -1, "xmin<0 rejected");
    CHECK(ei_write_table(file, 3.0, 2.0, 5) == -1, "inverted interval rejected");
    CHECK(ei_write_table(file, NAN, 2.0, 5) == -1, "NaN endpoint rejected");
    CHECK(ei_write_table(file, 1.0, INFINITY, 5) == -1,
          "infinite endpoint rejected");
    CHECK(ei_write_table(file, 1.0, 2.0, 0) == -1, "count=0 rejected");
    CHECK(ei_write_table(file, 1.0, 2.0, EI_MAX_POINTS + 1) == -1,
          "excessive count rejected");
    CHECK(ftell(file) == 0, "invalid arguments produce no output");
    fclose(file);
}

int main(void)
{
    test_reference_values();
    test_zero_and_domain();
    test_close_to_zero_of_ei();
    test_table_three_points();
    test_single_point_and_equal_endpoints();
    test_invalid_tables();

    printf("Test groups: %u; checks: %u; failures: %u\n", cases, checks, failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
