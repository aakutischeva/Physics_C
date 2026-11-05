#include "statistics.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_result(const char *name, double value)
{
    if (isnan(value))
        printf("%s: undefined\n", name);
    else
        printf("%s: %.17g\n", name, value);
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s data.txt\n", argv[0]);
        fprintf(stderr, "Use '-' instead of a filename to read standard input.\n");
        return EXIT_FAILURE;
    }

    int from_stdin = strcmp(argv[1], "-") == 0;
    FILE *input = from_stdin ? stdin : fopen(argv[1], "r");
    if (input == NULL) {
        perror("Cannot open input file");
        return EXIT_FAILURE;
    }

    Statistics stats;
    StatsStatus status = statistics_from_stream(input, &stats);
    if (!from_stdin && fclose(input) != 0 && status == STATS_OK)
        status = STATS_READ_ERROR;

    if (status != STATS_OK && status != STATS_SINGLE_VALUE
        && status != STATS_ZERO_VARIANCE) {
        fprintf(stderr, "Error: %s.\n", stats_status_message(status));
        return EXIT_FAILURE;
    }

    printf("Count: %u\n", stats.count);
    print_result("Mean", stats.mean);
    print_result("Variance (N-1)", stats.variance);
    print_result("Standard deviation", stats.stddev);
    print_result("Skewness", stats.skewness);
    print_result("Excess kurtosis", stats.excess_kurtosis);

    if (status != STATS_OK)
        fprintf(stderr, "Note: %s.\n", stats_status_message(status));

    return EXIT_SUCCESS;
}
