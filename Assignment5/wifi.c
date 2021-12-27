// PX425 2021 Assignment 5
// Determining the Percolation Threshold for an ad-hoc Wireless
// Network in a very large 3D model of a spherical space station
//
// Starter code for optimisation and parallelisation
// N. Hine, November 2021
// phuwcs, December 2021
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include "mt19937ar.h"

// Typedefs

typedef struct Router {
    /// Position from centre
    double x, y, z;
    /// Router radius
    double m, r;
    /// Cluster number
    int cluster;
} Router;

/// Domain decomposition of a volume filled with routers, containing sizes, number of cells,
/// pointers to the routers in each cell, total number of routers, total number of clusters,
/// and whether the clusters span the space.
typedef struct CellDomain {
    /// Size
    double S, Lx, Ly, Lz;
    /// Number of cells and sizes of cells
    int nx, ny, nz, nc;
    struct Router*** cell_rtr;
    /// Number of routers in each cell
    int* cell_nrtr;
    /// Total number of routers
    int nrtr;
    /// Cluster number
    int ncluster;
    /// True if clusters span the space
    bool spanning_cluster;
} CellDomain;

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

// Function Prototypes

void parse_int(int* addr, char* err_name);

void parse_double(double* addr, char* err_name);

int strtoi(const char* nptr, char** endptr, int base);

bool is_parsed_int_invalid(const char* startptr, const char* endptr, int value);

bool is_parsed_double_invalid(const char* startptr, const char* endptr, double value);

int read_args(int argc, char** argv, Args* args);

void find_cluster(int Nrtr, Router** Rtr, int i);

bool connected(Router* ra, Router* rb);

int find_spanning_cluster(CellDomain* dom);

bool merge_clusters(int nra, Router** ra, int nrb, Router** rb);

int count_clusters(CellDomain* dom);

void generate_routers(int* Nrtr, Router** Rtr, double S, double R, double P);

void create_domain_decomp(CellDomain* dom, int Nrtr, Router* Rtr, int nx, int ny, int nz);

void find_all_clusters(CellDomain* dom);

void destroy_domain_decomp(CellDomain* dom);

// Main Routine
int main(int argc, char** argv) {
    // Default Values
    Args args = {
            .router_radius = 1.0,
            .space_station_initial_size = 20.0,
            .space_station_size_increment = 5.0,
            .number_space_station_size_increments = 5,
            .volume_fraction_initial = 0.34,
            .volume_fraction_increment = 0.015,
            .number_volume_fraction_increments = 1,
            .use_current_time_as_seed = false,
    };

    // Read command line arguments
    int ra = read_args(argc, argv, &args);
    if (ra == 0) {
        printf("# Command line arguments successfully processed\n");
    } else {
        printf("Error Processing command line arguments (error code %i)\n", ra);
        return EXIT_FAILURE;
    }

    // Create bindings
    double R = args.router_radius;

    double Sinit = args.space_station_initial_size;
    double deltaS = args.space_station_size_increment;
    int nS = args.number_space_station_size_increments;

    double Pinit = args.volume_fraction_initial;
    double deltaP = args.volume_fraction_increment;
    int nP = args.number_volume_fraction_increments;

    bool seedtime = args.use_current_time_as_seed;

    int cellmin = 1, cellmax = 1;

    // Seed random number generator
    unsigned long seed = 20350292;
    if (seedtime) {
        seed = time(NULL);
        printf("# Using time-based random seed %ld\n", seed);
    }
    init_genrand(seed);

    // Router-related variables
    int Nrtr = 0;
    Router* Rtr = NULL;
    // Cell decomposition information
    CellDomain dom;
    // Number of runs for this invocation of the program
    int nruns = nP * nS;

    // Loop over values of filling fraction P and system size space_station_initial_size
    for (int irun = 0; irun < nruns; irun++) {
        // Find values of iP and iS for this run
        int iP = irun % nP;
        int iS = (irun - iP) / nP;
        // Find values of P and space_station_initial_size for this run
        double P = Pinit + iP * deltaP;
        double S = Sinit + iS * deltaS;
        dom.S = S;
        dom.Lx = 2.0 * S;
        dom.Ly = 2.0 * S;
        dom.Lz = 2.0 * S;
        // Generate randomly-placed routers in the domain
        generate_routers(&Nrtr, &Rtr, S, R, P);
        // Output sizes and volume fraction (+newline if multiple cell sizes)
        printf("space_station_initial_size = %6.2f P = %8.6f ", S, P);
        if (cellmin != cellmax) { printf("\n"); }
        // Loop over domain decomposition grid sizes
        for (int i = cellmin; i <= cellmax; i++) {
            if (cellmin != cellmax) { printf("ncells = %3d ", i); }
            // Initialise the domain decomposition structure
            create_domain_decomp(&dom, Nrtr, Rtr, i, i, i);
            // Find clusters in cells, merge between cells, count clusters
            // and find spanning cluster if it exists
            find_all_clusters(&dom);
            // Write results to stdout
            printf(": %6d clusters, ", dom.ncluster);
            if (dom.spanning_cluster == 0) {
                printf("none spanning\n");
            } else {
                printf("%7d spans\n", dom.spanning_cluster);
            }
            // remove storage associated with domain decomposition and
            // reset Router cluster values
            destroy_domain_decomp(&dom);
            for (int j = 0; j < Nrtr; j++) {
                Rtr[j].cluster = 0;
            }
        }
        free(Rtr);
    }
    return EXIT_SUCCESS;
}

/// Set up number of routers for a given volume fraction, and
/// generate random coordinates for each one, setting cluster
/// value to 0 to indicate cluster not found yet.
void generate_routers(int* Nrtr, Router** Rtr, double S, double R, double P) {
    *Nrtr = (int) (S * S * S * P / (R * R * R));
    *Rtr = calloc(*Nrtr, sizeof(Router));
    double mag, x, y, z;
    for (int i = 0; i < (*Nrtr); i++) {
        // Find random coordinates inside the sphere of radius space_station_initial_size
        // whose centre is at (space_station_initial_size,space_station_initial_size,space_station_initial_size)
        mag = S * 4.0;
        while (mag > S) {
            x = genrand() * S * 2.0;
            y = genrand() * S * 2.0;
            z = genrand() * S * 2.0;
            mag = sqrt((x - S) * (x - S) + (y - S) * (y - S) + (z - S) * (z - S));
        }
        (*Rtr)[i].x = x;
        (*Rtr)[i].y = y;
        (*Rtr)[i].z = z;
        (*Rtr)[i].m = mag;
        (*Rtr)[i].r = R;
        (*Rtr)[i].cluster = 0;
    }
}

/// Set up a cartesian grid of cells as a "domain decomposition" for helping
/// find clusters of routers which span the space. Assign each Router to
/// a cell, and set up lists of pointers to the routers in each cell.
void create_domain_decomp(CellDomain* dom, int Nrtr, Router* Rtr, int nx, int ny, int nz) {

    dom->nx = nx;
    dom->ny = ny;
    dom->nz = nz;
    dom->nc = nx * ny * nz;
    dom->cell_nrtr = calloc(dom->nx * dom->ny * dom->nz, sizeof(int));
    dom->cell_rtr = calloc(dom->nx * dom->ny * dom->nz, sizeof(Router*));
    dom->nrtr = Nrtr;
    dom->spanning_cluster = 0;

    // count the routers associated with each cell
    for (int i = 0; i < Nrtr; i++) {
        int ix = floor(dom->nx * Rtr[i].x / dom->Lx);
        int iy = floor(dom->ny * Rtr[i].y / dom->Ly);
        int iz = floor(dom->nz * Rtr[i].z / dom->Lz);
        int ic = ix + dom->ny * iy + dom->ny * dom->nz * iz;
        dom->cell_nrtr[ic]++;
    }

    // allocate memory for storage of each cell's set of pointers to routers,
    // then re-set the counts of numbers of routers to zero
    for (int ic = 0; ic < dom->nc; ic++) {
        dom->cell_rtr[ic] = calloc(dom->cell_nrtr[ic], sizeof(Router*));
        dom->cell_nrtr[ic] = 0;
    }

    // fill up the list of routers associated with each cell, counting them
    // again as we go for indexing purposes
    for (int i = 0; i < Nrtr; i++) {
        int ix = floor(dom->nx * Rtr[i].x / dom->Lx);
        int iy = floor(dom->ny * Rtr[i].y / dom->Ly);
        int iz = floor(dom->nz * Rtr[i].z / dom->Lz);
        int ic = ix + dom->ny * iy + dom->ny * dom->nz * iz;
        dom->cell_rtr[ic][dom->cell_nrtr[ic]] = &Rtr[i];
        dom->cell_nrtr[ic]++;
    }
}

/// Deallocate storage associated with a domain decomposition struct.
void destroy_domain_decomp(CellDomain* dom) {
    int ic;
    for (ic = 0; ic < dom->nc; ic++) {
        free(dom->cell_rtr[ic]);
    }
    free(dom->cell_rtr);
    free(dom->cell_nrtr);
}

/// Identify all clusters in a domain decomposition structure.
void find_all_clusters(CellDomain* dom) {

    clock_t t1 = clock();
    int cl = 1;
    // loop over all cells of the domain decomposition, then loop over the routers
    // in each cell. If cluster has not yet been identified, find all the
    // connected routers within this cell's list of routers (fast!)
    for (int ic = 0; ic < dom->nc; ic++) {
        for (int i = 0; i < dom->cell_nrtr[ic]; i++) {
            if (dom->cell_rtr[ic][i]->cluster == 0) {
                dom->cell_rtr[ic][i]->cluster = cl;
                find_cluster(dom->cell_nrtr[ic], dom->cell_rtr[ic], i);
                cl++;
            }
        }
    }
    clock_t t2 = clock();

    // merge clusters between cells if they are connected. Start from first
    // cell and move outwards, checking the "outward" half of the set of nearest
    // neighbour cells. Always retain lower numbered cluster to prevent circular
    // "flows" of cluster value
    int changed = 1;
    // keep repeating loop until nothing changes any more
    while (changed > 0) {
        changed = 0;
        for (int ic = 0; ic < dom->nx * dom->ny * dom->nz; ic++) {
            // loop over 13 of the 26 "nearest neighbour" cells on the cubic
            // lattice, ie the outward half - otherwise double counting will
            // occur and waste time
            int ix = ic % dom->nx;
            int iy = ((ic - ix) / dom->nx) % dom->ny;
            int iz = (ic - ix - dom->ny * iy) / (dom->nx * dom->ny);
            for (int neighb = 14; neighb < 27; neighb++) {
                // modulo arithmetic to find dx,dy,dz
                int dx = neighb % 3 - 1;
                int dy = ((neighb - dx) / 3) % 3 - 1;
                int dz = (neighb - dx - 3 * dy) / 9 - 1;
                // prevent checking beyond limits of cell grid
                if ((ix + dx >= dom->nx) || (iy + dy >= dom->ny) || (iz + dz >= dom->nz)) { continue; }
                if ((ix + dx < 0) || (iy + dy < 0) || (iz + dz < 0)) { continue; }
                // find index of neighbour cell
                int icp = ic + dx + dom->ny * dy + dom->ny * dom->nz * dz;
                changed = changed + merge_clusters(dom->cell_nrtr[ic], dom->cell_rtr[ic], dom->cell_nrtr[icp],
                                                   dom->cell_rtr[icp]);
            }
        }
    }
    clock_t t3 = clock();
    dom->ncluster = count_clusters(dom);
    dom->spanning_cluster = find_spanning_cluster(dom);
    clock_t t4 = clock();
    // print timings - you may need to disable this in parallel if
    // it is causing problems
    printf("(%8.4f %8.4f %8.4f sec) ", (double) (t2 - t1) / (double) CLOCKS_PER_SEC,
           (double) (t3 - t2) / (double) CLOCKS_PER_SEC, (double) (t4 - t3) / (double) CLOCKS_PER_SEC);
}

/// Find all the "connected" routers in a list, starting from a specific Router i,
/// and sets their cluster value to match that of the starting Router.
void find_cluster(int Nrtr, Router** Rtr, int i) {
    for (int j = 0; j < Nrtr; j++) {
        if (connected(Rtr[i], Rtr[j]) && (Rtr[j]->cluster != Rtr[i]->cluster)) {
            Rtr[j]->cluster = Rtr[i]->cluster;
            find_cluster(Nrtr, Rtr, j);
        }
    }
}

/// Check if two routers are "connected", ie their separation is less than the sum of their radii.
bool connected(Router* ra, Router* rb) {
    double separation =
            (ra->x - rb->x) * (ra->x - rb->x) + (ra->y - rb->y) * (ra->y - rb->y) + (ra->z - rb->z) * (ra->z - rb->z);
    double sum_of_radii = (ra->r + rb->r) * (ra->r + rb->r);
    if (separation <= sum_of_radii) { return true; }
    else { return false; }
}

/// Merge the clusters associated with two lists of routers if pairs of them are closer than the sum of their radii.
///
/// Retains lower-numbered cluster and converts higher-numbered cluster to match.
bool merge_clusters(int nra, Router** ra, int nrb, Router** rb) {
    bool changed = false;
    // Loop over routers i, j in the two cells
    for (int i = 0; i < nra; i++) {
        for (int j = 0; j < nrb; j++) {
            // search for pairs of routers whose cluster values are not equal
            if (ra[i]->cluster != rb[j]->cluster) {
                // if they are closer than the sum of their radii ..
                if (connected(ra[i], rb[j])) {
                    // convert whichever cluster has the higher index to match the
                    // lower index
                    if (rb[j]->cluster > ra[i]->cluster) {
                        // cluster in cell A has lower index
                        int cl = rb[j]->cluster;
                        for (int k = 0; k < nrb; k++) {
                            if (rb[k]->cluster == cl) {
                                rb[k]->cluster = ra[i]->cluster;
                            }
                        }
                    }  // else cluster in cell B has lower index
                    else {
                        int cl = ra[i]->cluster;
                        for (int k = 0; k < nra; k++) {
                            if (ra[k]->cluster == cl) {
                                ra[k]->cluster = rb[j]->cluster;
                            }
                        }
                    }
                    // set flag to remember that something changed within this call so
                    // that we can halt merging once nothing changes any more
                    changed = true;
                }
            }
        }
    }
    return changed;
}

///Count the number of unique clusters in the list of routers
int count_clusters(CellDomain* dom) {
    int count = 0;
    int* counted = calloc(dom->nrtr, sizeof(int));
    for (int ic = 0; ic < dom->nc; ic++) {
        for (int i = 0; i < dom->cell_nrtr[ic]; i++) {
            bool found = false;
            // check if we have already counted this cluster
            for (int j = 0; j < count; j++) {
                if (dom->cell_rtr[ic][i]->cluster == counted[j]) {
                    found = true;
                    break;
                }
            }
            // if not, add it to the list and increment the count
            if (!found) {
                counted[count] = dom->cell_rtr[ic][i]->cluster;
                count++;
            }
        }
    }
    free(counted);
    return count;
}

///check if there are clusters extending to the surface in each octant
int find_spanning_cluster(CellDomain* dom) {
    int* octant_check_list[8];

    int spanning_cluster = 0;

    // Set up and initialise storage for finding spanning cluster
    for (int ioct = 0; ioct < 8; ioct++) {
        octant_check_list[ioct] = calloc(dom->ncluster, sizeof(int));
        for (int i = 0; i < dom->ncluster; i++) {
            octant_check_list[ioct][i] = -1;
        }
    }

    //Loop over cells
    for (int ic = 0; ic < dom->nc; ic++) {
        // Loop over routers in each cell
        for (int i = 0; i < dom->cell_nrtr[ic]; i++) {
            double x = dom->cell_rtr[ic][i]->x;
            double y = dom->cell_rtr[ic][i]->y;
            double z = dom->cell_rtr[ic][i]->z;
            double m = dom->cell_rtr[ic][i]->m;
            double r = dom->cell_rtr[ic][i]->r;
            int cl = dom->cell_rtr[ic][i]->cluster;
            // If Router touches the sphere edge...
            if (m + r > dom->S) {
                int ioct = 0;
                // Calculate which octant Router is in
                if (x > dom->S) { ioct++; }
                if (y > dom->S) { ioct += 2; }
                if (z > dom->S) { ioct += 4; }
                // Check if this cluster has been found for this octant
                for (int icl = 0; icl < dom->ncluster; icl++) {
                    // If we reach the end of the list of clusters, add this one at
                    // the end and stop looking
                    if (octant_check_list[ioct][icl] == -1) {
                        octant_check_list[ioct][icl] = cl;
                        break;
                    }
                    // If this cluster has already been listed, stop looking
                    if (octant_check_list[ioct][icl] == cl) { break; }
                }
            }
        }
    }
    // Check if the same cluster appears in all 8 lists
    for (int icl = 0; icl < dom->ncluster; icl++) {
        int sum = 0;
        // Get next cluster from octant 0, quit if blank
        int cl = octant_check_list[0][icl];
        if (cl == -1) { break; }
        // Check all 8 octants to check this cluster appears in each
        for (int ioct = 0; ioct < 8; ioct++) {
            for (int jcl = 0; jcl < dom->ncluster; jcl++) {
                if (octant_check_list[ioct][jcl] == cl) {
                    sum++;
                    break;
                }
            }
        }
        if (sum == 8) { spanning_cluster++; }
    }
    // Free storage associated with spanning cluster check
    for (int ioct = 0; ioct < 8; ioct++) {
        free(octant_check_list[ioct]);
    }
    return spanning_cluster;
}

/// Parse Command line arguments
int read_args(int argc, char** argv, Args* args) {
    int c;
    opterr = 0;

    // Process all flags found in argv
    while ((c = getopt(argc, argv, "t:P:S:R:D:Q:N:M")) != -1) {
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
