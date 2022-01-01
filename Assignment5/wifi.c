// PX425 2021 Assignment 5
// Determining the Percolation Threshold for an ad-hoc Wireless
// Network in a very large 3D model of a spherical space station
//
// Starter code for optimisation and parallelisation
// N. Hine, November 2021
// phuwcs, December 2021
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include "mt19937ar.h"
#include "input.h"
#include "output.h"

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

// Function Prototypes

void find_cluster(int Nrtr, Router** Rtr, int i);

bool connected(Router* ra, Router* rb);

int find_spanning_cluster(CellDomain* dom);

bool merge_clusters(int nra, Router** ra, int nrb, Router** rb);

int count_clusters(CellDomain* dom);

void generate_routers(int* Nrtr, Router** Rtr, double S, double R, double P);

void create_domain_decomp(CellDomain* dom, int Nrtr, Router* Rtr, int nx, int ny, int nz);

void find_all_clusters(CellDomain* dom, DecompResults* decomp_results);

void destroy_domain_decomp(CellDomain* dom);

// Main Routine
int main(int argc, char** argv) {
    double start_time = omp_get_wtime();
    // Default Values
    Args args = {
            .router_radius                        = 1.0,
            .space_station_initial_size           = 20.0,
            .space_station_size_increment         = 5.0,
            .number_space_station_size_increments = 5,
            .volume_fraction_initial              = 0.34,
            .volume_fraction_increment            = 0.015,
            .number_volume_fraction_increments    = 1,
            .use_current_time_as_seed             = false,
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

    int cellmin = 20, cellmax = 20;

    // Seed random number generator
    unsigned long seed = 20350292;
    if (seedtime) {
        seed = time(NULL);
        printf("# Using time-based random seed %ld\n", seed);
    }
    init_genrand(seed);

    // Number of runs for this invocation of the program
    int nruns = nP * nS;

    Results results = {
            .nruns = nruns,
    };
    create_Results(&results);

    // Loop over values of filling fraction P and system size S
    for (int irun = 0; irun < nruns; irun++) {
        // Find values of iP and iS for this run
        int iP = irun % nP;
        int iS = (irun - iP) / nP;
        // Find values of P and S for this run
        double P = Pinit + iP * deltaP;
        double S = Sinit + iS * deltaS;
        // Cell decomposition information
        CellDomain dom = {
                .S = S,
                .Lx = 2.0 * S,
                .Ly = 2.0 * S,
                .Lz = 2.0 * S,
        };
        // Router-related variables
        int Nrtr;
        // Generate randomly-placed routers in the domain
        Router* Rtr;
        generate_routers(&Nrtr, &Rtr, S, R, P);
        // Output sizes and volume fraction
        RunResults* run_results = &results.run_results[irun];
        run_results->S = S;
        run_results->P = P;
        run_results->n_cell_sizes = cellmax - cellmin + 1;
        create_RunResults(run_results);
        // Loop over domain decomposition grid sizes
        for (int i = cellmin; i <= cellmax; i++) {
            DecompResults* decomp_results = &run_results->decomp_results[i - cellmin];
            decomp_results->ncells = i;
            // Initialise the domain decomposition structure
            create_domain_decomp(&dom, Nrtr, Rtr, i, i, i);
            // Find clusters in cells, merge between cells, count clusters
            // and find spanning cluster if it exists
            find_all_clusters(&dom, decomp_results);
            // Save results
            decomp_results->n_clusters = dom.ncluster;
            decomp_results->spanning_cluster = dom.spanning_cluster;
            // remove storage associated with domain decomposition and
            // reset Router cluster values
            destroy_domain_decomp(&dom);
            for (int j = 0; j < Nrtr; j++) {
                Rtr[j].cluster = 0;
            }
        }
        free(Rtr);
    }

    print_Results(&results);
    destroy_Results(&results);

    double end_time = omp_get_wtime();
    printf("\nTotal time: %6.4f\n", end_time - start_time);

    return EXIT_SUCCESS;
}

/// Set up number of routers for a given volume fraction, and
/// generate random coordinates for each one, setting cluster
/// value to 0 to indicate cluster not found yet.
void generate_routers(int* Nrtr, Router** Rtr, double S, double R, double P) {
    *Nrtr = (int) (S * S * S * P / (R * R * R));
    *Rtr = calloc(*Nrtr, sizeof(Router));
    if (*Rtr == NULL) {
        fprintf(stderr, "Failed to allocate %i elements of size %lu\n", *Nrtr, sizeof(Router));
        abort();
    }
    for (int i = 0; i < (*Nrtr); i++) {
        double mag = S * 4.0;
        double x = 0.0, y = 0.0, z = 0.0;
        // Find random coordinates inside the sphere of radius S
        // whose centre is at (S, S, S)
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
    dom->cell_nrtr = calloc(dom->nc, sizeof(int));
    dom->cell_rtr = calloc(dom->nc, sizeof(Router*));
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
void find_all_clusters(CellDomain* dom, DecompResults* decomp_results) {

    double t1 = omp_get_wtime();
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

    double t2 = omp_get_wtime();

    // merge clusters between cells if they are connected. Start from first
    // cell and move outwards, checking the "outward" half of the set of nearest
    // neighbour cells. Always retain lower numbered cluster to prevent circular
    // "flows" of cluster value
    bool changed = true;
    int dxs[27], dys[27], dzs[27];
    for (int i = 14; i < 27; i++) {
        dxs[i] = i % 3 - 1;
        dys[i] = ((i - dxs[i]) / 3) % 3 - 1;
        dzs[i] = (i - dxs[i] - 3 * dys[i]) / 9 - 1;
    }

    // keep repeating loop until nothing changes any more
    while (changed) {
        changed = false;
#pragma omp parallel for default(none) shared(dom, dxs, dys, dzs) reduction(||: changed)
        for (int ic = 0; ic < dom->nx * dom->ny * dom->nz; ic++) {
            // loop over 13 of the 26 "nearest neighbour" cells on the cubic
            // lattice, ie the outward half - otherwise double counting will
            // occur and waste time
            int ix = ic % dom->nx;
            int iy = ((ic - ix) / dom->nx) % dom->ny;
            int iz = (ic - ix - dom->ny * iy) / (dom->nx * dom->ny);
            for (int neighb = 14; neighb < 27; neighb++) {
                // modulo arithmetic to find dx, dy, dz
                int dx = dxs[neighb];
                int dy = dys[neighb];
                int dz = dzs[neighb];
                // prevent checking beyond limits of cell grid
                if ((ix + dx >= dom->nx) || (iy + dy >= dom->ny) || (iz + dz >= dom->nz)) { continue; }
                if ((ix + dx < 0) || (iy + dy < 0) || (iz + dz < 0)) { continue; }
                // find index of neighbour cell
                int icp = ic + dx + dom->ny * dy + dom->ny * dom->nz * dz;
                changed |= merge_clusters(dom->cell_nrtr[ic], dom->cell_rtr[ic], dom->cell_nrtr[icp],
                                          dom->cell_rtr[icp]);
            }
        }
    }
    double t3 = omp_get_wtime();
    dom->ncluster = count_clusters(dom);
    dom->spanning_cluster = find_spanning_cluster(dom);
    double t4 = omp_get_wtime();
    DecompResults_add_times(decomp_results, t1, t2, t3, t4);
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
            if (ra[i]->cluster == rb[j]->cluster) { continue; }
            // if they are closer than the sum of their radii ..
            if (!connected(ra[i], rb[j])) { continue; }
            // convert whichever cluster has the higher index to match the lower index
            Router** higher, ** lower;
            int high_index, low_index, n_high;
            if (rb[j]->cluster > ra[i]->cluster) {
                // cluster in cell A has lower high_index
                lower = ra, higher = rb;
                low_index = i, high_index = j;
                n_high = nrb;
            } else {
                // else cluster in cell B has lower high_index
                lower = rb, higher = ra;
                low_index = j, high_index = i;
                n_high = nra;
            }
            int cl = higher[high_index]->cluster;
            for (int k = 0; k < n_high; k++) {
                if (higher[k]->cluster != cl) { continue; }
                higher[k]->cluster = lower[low_index]->cluster;
            }
            // set flag to remember that something changed within this call so
            // that we can halt merging once nothing changes any more
            changed = true;
        }
    }
    return changed;
}

/// Count the number of unique clusters in the list of routers
int count_clusters(CellDomain* dom) {
    // The max number of clusters is the number of routers (each router in its own cluster)
    bool* counted = calloc(dom->nrtr, sizeof(bool));
    // This is a relaxed loop since it's just setting things to true
#pragma omp parallel for default(none) shared(dom, counted)
    for (int ic = 0; ic < dom->nc; ic++) {
        for (int i = 0; i < dom->cell_nrtr[ic]; i++) {
            // Mark this cluster as having been found
            int cluster = dom->cell_rtr[ic][i]->cluster;
            counted[cluster] = true;
        }
    }

    // Count the number of clusters that have been marked as found
    int count = 0;
    for (int i = 0; i < dom->nrtr; i++) {
        count += counted[i];
    }
    free(counted);
    return count;
}

/// Check if there are clusters extending to the surface in each octant
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
#pragma omp parallel for default(none) shared(dom, octant_check_list)
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
#pragma omp parallel for default(none) shared(dom, octant_check_list, cl) reduction(+: sum)
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

