#include "statistics.h"

#include <float.h>
#include <limits.h>
#include <math.h>

/* Streaming updates of central moments M2, M3, M4 (Pébay/Welford).
   Mi = sum_j (x_j - mean)^i. No rewind or data array is needed. */
StatsStatus statistics_from_stream(FILE *input, Statistics *result)
{
    if (input == NULL || result == NULL)
        return STATS_INVALID_ARGUMENT;

    *result = (Statistics){0, NAN, NAN, NAN, NAN, NAN};

    unsigned int count = 0;
    /* Subtract the first observation before updating the moments.  A
       large common offset (e.g. 1e12) then cannot degrade the mean update. */
    long double origin = 0.0L;
    long double mean = 0.0L;
    long double m2 = 0.0L;
    long double m3 = 0.0L;
    long double m4 = 0.0L;

    for (;;) {
        double value;
        int scanned = fscanf(input, "%lf", &value);
        if (scanned == EOF) {
            if (ferror(input))
                return STATS_READ_ERROR;
            break;
        }
        if (scanned != 1 || !isfinite(value))
            return STATS_INVALID_DATA;
        if (count == UINT_MAX)
            return STATS_COUNT_OVERFLOW;

        if (count == 0)
            origin = (long double)value;
        long double centered = (long double)value - origin;
        long double previous_count = (long double)count;
        long double n = (long double)++count;
        long double delta = centered - mean;
        long double delta_n = delta / n;
        long double delta_n2 = delta_n * delta_n;
        long double term1 = delta * delta_n * previous_count;

        /* Use the previous moments in all three updates. */
        m4 += term1 * delta_n2 * (n * n - 3.0L * n + 3.0L)
              + 6.0L * delta_n2 * m2 - 4.0L * delta_n * m3;
        m3 += term1 * delta_n * (n - 2.0L) - 3.0L * delta_n * m2;
        m2 += term1;
        mean += delta_n;

        if (!isfinite(mean) || !isfinite(m2)
            || !isfinite(m3) || !isfinite(m4))
            return STATS_NUMERIC_RANGE;
    }

    result->count = count;
    if (count == 0)
        return STATS_EMPTY;

    result->mean = (double)(origin + mean);
    if (!isfinite(result->mean))
        return STATS_NUMERIC_RANGE;

    if (count == 1)
        return STATS_SINGLE_VALUE;

    long double variance = m2 / (long double)(count - 1);
    result->variance = (double)variance;
    if (!isfinite(result->variance)
        || (variance > 0.0L && result->variance == 0.0))
        return STATS_NUMERIC_RANGE;

    if (m2 == 0.0L) {
        result->stddev = 0.0;
        return STATS_ZERO_VARIANCE;
    }

    long double sigma = sqrtl(variance);
    long double skewness = (m3 / (long double)count)
                           / (variance * sigma);
    long double excess = (m4 / (long double)count)
                         / (variance * variance) - 3.0L;

    result->stddev = (double)sigma;
    result->skewness = (double)skewness;
    result->excess_kurtosis = (double)excess;

    if (!isfinite(result->stddev) || result->stddev == 0.0
        || !isfinite(result->skewness)
        || !isfinite(result->excess_kurtosis))
        return STATS_NUMERIC_RANGE;

    return STATS_OK;
}

const char *stats_status_message(StatsStatus status)
{
    switch (status) {
    case STATS_OK: return "all statistics are defined";
    case STATS_EMPTY: return "the file contains no numbers";
    case STATS_SINGLE_VALUE: return "N = 1: sample variance is undefined";
    case STATS_ZERO_VARIANCE: return "all values coincide: skewness and kurtosis are undefined";
    case STATS_INVALID_DATA: return "invalid or non-finite number in input";
    case STATS_READ_ERROR: return "file read error";
    case STATS_COUNT_OVERFLOW: return "too many observations";
    case STATS_NUMERIC_RANGE: return "result outside the representable double range";
    case STATS_INVALID_ARGUMENT: return "invalid function argument";
    }
    return "unknown error";
}
