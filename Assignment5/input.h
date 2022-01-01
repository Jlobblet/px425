#ifndef PX425_INPUT_H
#define PX425_INPUT_H

#include <stdbool.h>

typedef struct Args {
    /// Radius of each Router
    double router_radius;
    /// Space station size starting value
    double space_station_initial_size;
    /// Increment of space station size
    double space_station_size_increment;
    /// Number of increments of space station size
    int number_space_station_size_increments;
    /// Volume fraction starting value
    double volume_fraction_initial;
    /// Increment of volume fraction
    double volume_fraction_increment;
    /// Number of increments of volume fraction
    int number_volume_fraction_increments;
    /// If true, use current time as random seed
    bool use_current_time_as_seed;
} Args;

int read_args(int argc, char** argv, Args* args);

#endif //PX425_INPUT_H
