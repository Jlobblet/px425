#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <errno.h>
#include <ctype.h>
#include "input.h"

void parse_int(int* addr, char* err_name);

void parse_double(double* addr, char* err_name);

int strtoi(const char* nptr, char** endptr, int base);

bool is_parsed_int_invalid(const char* startptr, const char* endptr, int value);

bool is_parsed_double_invalid(const char* startptr, const char* endptr, double value);

/// Parse Command line arguments
int read_args(int argc, char** argv, Args* args) {
    int c;
    opterr = 0;

    // Process all flags found in argv
    while ((c = getopt(argc, argv, "tP:S:R:D:Q:N:M:")) != -1) {
        switch (c) {
            case 't':
                args->use_current_time_as_seed = true;
                break;
            case 'P':
                parse_double(&args->volume_fraction_initial, "Volume fraction");
                break;
            case 'S':
                parse_double(&args->space_station_initial_size, "Station size");
                break;
            case 'R':
                parse_double(&args->router_radius, "Router radius");
                break;
            case 'D':
                parse_double(&args->space_station_size_increment, "Station step size");
                break;
            case 'Q':
                parse_double(&args->volume_fraction_increment, "Volume fraction step");
                break;
            case 'N':
                parse_int(&args->number_space_station_size_increments, "Number of station size steps");
                break;
            case 'M':
                parse_int(&args->number_volume_fraction_increments, "Number of volume fraction steps");
                break;
            case '?':
                if ((optopt == 'R') || (optopt == 'P') || (optopt == 'S')) {
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                } else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                } else {
                    fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
                }
                return 1;
            default:
                abort();
        }
    }

    // List unrecognised arguments
    for (int index = optind; index < argc; index++) {
        printf("Non-option argument %s\n", argv[index]);
    }

    return opterr;
}

void parse_int(int* addr, char* err_name) {
    char* endptr;
    *addr = strtoi(optarg, &endptr, 10);
    if (is_parsed_int_invalid(optarg, endptr, *addr)) {
        opterr = 1;
        fprintf(stderr, "%s argument could not be read: %s\n", err_name, optarg);
    }
}

void parse_double(double* addr, char* err_name) {
    char* endptr;
    *addr = strtod(optarg, &endptr);
    if (is_parsed_double_invalid(optarg, endptr, *addr)) {
        opterr = 1;
        fprintf(stderr, "%s argument could not be read: %s\n", err_name, optarg);
    }
}

int strtoi(const char* nptr, char** endptr, int base) {
    long parsed = strtol(nptr, endptr, base);
    if (parsed > INT_MAX) {
        errno = ERANGE;
        return INT_MAX;
    }
    if (parsed < INT_MIN) {
        errno = ERANGE;
        return INT_MIN;
    }
    return (int) parsed;
}

bool is_parsed_int_invalid(const char* startptr, const char* endptr, int value) {
    return startptr == endptr
           || (value == INT_MAX && errno == ERANGE)
           || value <= 0;
}

bool is_parsed_double_invalid(const char* startptr, const char* endptr, double value) {
    return startptr == endptr
           || isinf(value)
           || isnan(value)
           || (value == HUGE_VAL && errno == ERANGE)
           || value <= 0.0;
}
