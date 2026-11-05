#include "ei.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define INPUT_SIZE 256

/* Read a complete line, rejecting embedded NULs and overly long lines.
   Discard the rest of a rejected line so it cannot become the next answer. */
static int read_line(char *buffer, size_t capacity)
{
    size_t used = 0;
    int ch;
    int invalid = 0;

    while ((ch = getchar()) != EOF && ch != '\n') {
        if (ch == '\0' || used + 1 >= capacity)
            invalid = 1;
        if (!invalid)
            buffer[used++] = (char)ch;
    }

    buffer[used] = '\0';
    return !invalid && (ch == '\n' || used != 0);
}

static void skip_spaces(const char **cursor)
{
    while (isspace((unsigned char)**cursor))
        ++*cursor;
}

static int end_of_input(const char *cursor)
{
    skip_spaces(&cursor);
    return *cursor == '\0';
}

static int parse_double(const char **cursor, double *value)
{
    skip_spaces(cursor);
    char *end;
    double parsed = strtod(*cursor, &end);
    if (end == *cursor || !isfinite(parsed))
        return 0;
    *value = parsed;
    *cursor = end;
    return 1;
}

/* The limit is 1,000,000, so unsigned int and the portable %u are enough.
   Parse the ENTIRE integer token: scanf("%u") accepts 3.5 and may wrap. */
static int parse_count(const char **cursor, unsigned int *value)
{
    skip_spaces(cursor);
    if (!isdigit((unsigned char)**cursor))
        return 0; /* Reject signs and non-integer input. */

    errno = 0;
    char *end;
    unsigned long parsed = strtoul(*cursor, &end, 10);
    if (errno != 0 || end == *cursor || parsed == 0
        || parsed > EI_MAX_POINTS)
        return 0;

    *value = (unsigned int)parsed;
    *cursor = end;
    return 1;
}

int main(void)
{
    unsigned int choice;
    char line[INPUT_SIZE];
    const char *cursor;

    printf("Section 2.2: exponential integral Ei(x), x > 0\n");
    printf("1 - evaluate Ei(x)\n");
    printf("2 - save a table to ei.dat\n");
    printf("Choice: ");
    fflush(stdout);

    if (!read_line(line, sizeof line)) {
        fprintf(stderr, "Error: enter 1 or 2.\n");
        return EXIT_FAILURE;
    }
    cursor = line;
    if (!parse_count(&cursor, &choice) || !end_of_input(cursor)
        || (choice != 1 && choice != 2)) {
        fprintf(stderr, "Error: enter 1 or 2.\n");
        return EXIT_FAILURE;
    }

    if (choice == 1) {
        double x;
        printf("Enter x (x > 0): ");
        fflush(stdout);
        if (!read_line(line, sizeof line)) {
            fprintf(stderr, "Error: x must be positive and finite.\n");
            return EXIT_FAILURE;
        }
        cursor = line;
        if (!parse_double(&cursor, &x) || !end_of_input(cursor) || x <= 0.0) {
            fprintf(stderr, "Error: x must be positive and finite.\n");
            return EXIT_FAILURE;
        }

        printf("Ei(%.17g) = %.17g\n", x, Ei(x));
        return EXIT_SUCCESS;
    }

    if (choice == 2) {
        double xmin, xmax;
        unsigned int count;
        printf("Enter xmin xmax number_of_points: ");
        fflush(stdout);
        if (!read_line(line, sizeof line)) {
            fprintf(stderr, "Error: require 0 < xmin <= xmax and 1..%u points.\n",
                    EI_MAX_POINTS);
            return EXIT_FAILURE;
        }
        cursor = line;
        if (!parse_double(&cursor, &xmin) || !parse_double(&cursor, &xmax)
            || !parse_count(&cursor, &count) || !end_of_input(cursor)
            || xmin <= 0.0 || xmax < xmin) {
            fprintf(stderr, "Error: require 0 < xmin <= xmax and 1..%u points.\n",
                    EI_MAX_POINTS);
            return EXIT_FAILURE;
        }

        FILE *out = fopen("ei.dat", "w");
        if (out == NULL) {
            perror("Error opening ei.dat");
            return EXIT_FAILURE;
        }

        int write_status = ei_write_table(out, xmin, xmax, count);
        int close_status = fclose(out);
        if (write_status != 0 || close_status != 0) {
            fprintf(stderr, "Error writing ei.dat.\n");
            return EXIT_FAILURE;
        }

        printf("Saved %u points to ei.dat (columns: x Ei(x)).\n", count);
        return EXIT_SUCCESS;
    }

    fprintf(stderr, "Error: unknown menu item.\n");
    return EXIT_FAILURE;
}
