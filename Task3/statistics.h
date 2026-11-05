#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdio.h>

typedef enum {
    STATS_OK = 0,
    STATS_EMPTY,
    STATS_SINGLE_VALUE,
    STATS_ZERO_VARIANCE,
    STATS_INVALID_DATA,
    STATS_READ_ERROR,
    STATS_COUNT_OVERFLOW,
    STATS_NUMERIC_RANGE,
    STATS_INVALID_ARGUMENT
} StatsStatus;

typedef struct {
    unsigned int count;
    double mean;
    double variance;
    double stddev;
    double skewness;
    double excess_kurtosis;
} Statistics;

/* One sequential pass, O(1) memory, no array of measurements.
   Variance uses N-1; skewness and excess kurtosis use N and that stddev.
   Undefined quantities are NAN for valid but degenerate datasets. */
StatsStatus statistics_from_stream(FILE *input, Statistics *result);

const char *stats_status_message(StatsStatus status);

#endif
