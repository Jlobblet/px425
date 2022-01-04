// PX425 2021 Assignment 5
// Determining the Percolation Threshold for an ad-hoc Wireless
// Network in a very large 3D model of a spherical space station
//
// Starter code for optimisation and parallelisation
// N. Hine, November 2021
// phuwcs, December 2021
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
#include "types.h"
#include "comms.h"

// Function Prototypes

void find_cluster(int Nrtr, Router** Rtr, int i);

bool connected(Router* ra, Router* rb);

int find_spanning_cluster(CellDomain* dom);

bool merge_clusters(int nra, Router** ra, int nrb, Router** rb);

int count_clusters(CellDomain* dom);

void find_all_clusters(CellDomain* dom, DecompResults* decomp_results);

void do_run(const Args* args, const MpiInfo* info, MPI_Request* request, int cellmin, int cellmax, RunResults* run_results, int irun);

// Main Routine
int main(int argc, char** argv) {
    MpiInfo info;
    comms_initialise(&argc, &argv, &info);
    const double start_time = omp_get_wtime();

    Args args = comms_read_input(argc, argv, &info);

    // Create bindings
    const int nS = args.number_space_station_size_increments;
    const int nP = args.number_volume_fraction_increments;

    bool seedtime = args.use_current_time_as_seed;

    const int cellmin = 20, cellmax = 20;

    // Seed random number generator
    unsigned long seed = 20350292;
    if (seedtime) {
        seed = time(NULL);
        printf("# Using time-based random seed %ld\n", seed);
    }
    init_genrand(seed);

    // Number of runs for this invocation of the program
    const int nruns = nP * nS;

    Results results = {
            .n_runs = nruns,
    };
    create_Results(&results);

    MPI_Request* requests = calloc(nruns, sizeof(MPI_Request));
    if (requests == NULL) { abort(); }
    MPI_RUNRESULTS = create_RunResults_datatype(cellmax - cellmin + 1);

    // Loop over values of filling fraction P and system size S
    for (int irun = 0; irun < nruns; irun++) {
        do_run(&args, &info, &requests[irun], cellmin, cellmax, &results.run_results[irun], irun);
    }

    if (info.my_rank == 0) {
        MPI_Status* statuses = calloc(nruns, sizeof(MPI_Status));
        if (statuses == NULL) { abort(); }
        MPI_Waitall(nruns, requests, statuses);
        print_Results(&results);
        const double end_time = omp_get_wtime();
        printf("\nTotal time: %6.4f\n", end_time - start_time);
        free(statuses);
    }

    free(requests);
    destroy_Results(&results);
    comms_finalise();

    return EXIT_SUCCESS;
}

void do_run(const Args* args, const MpiInfo* info, MPI_Request* request, int cellmin, int cellmax, RunResults* run_results, int irun) {
    // If this rank needs to do some work
    if (info->my_rank != 0 && (irun % info->n_processors) != info->my_rank) { return; }

    if (info->my_rank == 0) {
        MPI_Irecv(run_results, 1, MPI_RUNRESULTS, irun % info->n_processors, irun, MPI_COMM_WORLD, request);
    }
    if ((irun % info->n_processors) == info->my_rank) {
        // Find values of iP and iS for this run
        const int nP = args->number_volume_fraction_increments;
        const double Pinit = args->volume_fraction_initial;
        const double deltaP = args->volume_fraction_increment;
        const double Sinit = args->space_station_initial_size;
        const double deltaS = args->space_station_size_increment;
        const int iP = irun % nP;
        const int iS = (irun - iP) / nP;
        // Find values of P and S for this run
        const double P = Pinit + iP * deltaP;
        const double S = Sinit + iS * deltaS;
        // Output sizes and volume fraction
        run_results->S = S;
        run_results->P = P;
        run_results->n_cell_sizes = cellmax - cellmin + 1;
        create_RunResults(run_results);

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
        const double R = args->router_radius;
        generate_routers(&Nrtr, &Rtr, S, R, P);
        // Loop over domain decomposition grid sizes
//#pragma omp parallel default(none) shared(cellmin, cellmax, Nrtr, Rtr, run_results) private(dom)
        for (int i = cellmin; i <= cellmax; i++) {
            Router* Rtr_local = calloc(Nrtr, sizeof(Router));
            memcpy(Rtr_local, Rtr, Nrtr * sizeof(Router));
            DecompResults* decomp_results = &run_results->decomp_results[i - cellmin];
            decomp_results->ncells = i;
            // Initialise the domain decomposition structure
            create_domain_decomp(&dom, Nrtr, Rtr_local, i, i, i);
            // Find clusters in cells, merge between cells, count clusters
            // and find spanning cluster if it exists
            find_all_clusters(&dom, decomp_results);
            // Save results
            decomp_results->n_clusters = dom.n_clusters;
            decomp_results->spanning_cluster = dom.spanning_cluster;
            // remove storage associated with domain decomposition and
            // reset Router cluster values
            destroy_domain_decomp(&dom);
            free(Rtr_local);
        }
        MPI_Isend(run_results, 1, MPI_RUNRESULTS, 0, irun, MPI_COMM_WORLD, request);
        free(Rtr);
    }
}

/// Identify all clusters in a domain decomposition structure.
void find_all_clusters(CellDomain* dom, DecompResults* decomp_results) {
    double t1 = omp_get_wtime();
    int cl = 1;
    // loop over all cells of the domain decomposition, then loop over the routers
    // in each cell. If cluster has not yet been identified, find all the
    // connected routers within this cell's list of routers (fast!)
    for (int ic = 0; ic < dom->nc; ic++) {
        for (int i = 0; i < dom->cell_n_routers[ic]; i++) {
            if (dom->cell_routers[ic][i]->cluster == 0) {
                dom->cell_routers[ic][i]->cluster = cl;
                find_cluster(dom->cell_n_routers[ic], dom->cell_routers[ic], i);
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
                changed |= merge_clusters(dom->cell_n_routers[ic], dom->cell_routers[ic], dom->cell_n_routers[icp],
                                          dom->cell_routers[icp]);
            }
        }
    }
    double t3 = omp_get_wtime();
    dom->n_clusters = count_clusters(dom);
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
            (ra->x - rb->x) * (ra->x - rb->x)
            + (ra->y - rb->y) * (ra->y - rb->y)
            + (ra->z - rb->z) * (ra->z - rb->z);
    double sum_of_radii = (ra->r + rb->r) * (ra->r + rb->r);
    return separation <= sum_of_radii;
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
    bool* counted = calloc(dom->n_routers, sizeof(bool));
    // This is a relaxed loop since it's just setting things to true
#pragma omp parallel for default(none) shared(dom, counted)
    for (int ic = 0; ic < dom->nc; ic++) {
        for (int i = 0; i < dom->cell_n_routers[ic]; i++) {
            // Mark this cluster as having been found
            int cluster = dom->cell_routers[ic][i]->cluster;
            counted[cluster] = true;
        }
    }

    // Count the number of clusters that have been marked as found
    int count = 0;
    for (int i = 0; i < dom->n_routers; i++) {
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
        octant_check_list[ioct] = calloc(dom->n_clusters, sizeof(int));
        for (int i = 0; i < dom->n_clusters; i++) {
            octant_check_list[ioct][i] = -1;
        }
    }

    //Loop over cells
#pragma omp parallel for default(none) shared(dom, octant_check_list)
    for (int ic = 0; ic < dom->nc; ic++) {
        // Loop over routers in each cell
        for (int i = 0; i < dom->cell_n_routers[ic]; i++) {
            double x = dom->cell_routers[ic][i]->x;
            double y = dom->cell_routers[ic][i]->y;
            double z = dom->cell_routers[ic][i]->z;
            double m = dom->cell_routers[ic][i]->m;
            double r = dom->cell_routers[ic][i]->r;
            int cl = dom->cell_routers[ic][i]->cluster;
            // If Router touches the sphere edge...
            if (m + r > dom->S) {
                int ioct = 0;
                // Calculate which octant Router is in
                if (x > dom->S) { ioct++; }
                if (y > dom->S) { ioct += 2; }
                if (z > dom->S) { ioct += 4; }
                // Check if this cluster has been found for this octant
                for (int icl = 0; icl < dom->n_clusters; icl++) {
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
    for (int icl = 0; icl < dom->n_clusters; icl++) {
        int sum = 0;
        // Get next cluster from octant 0, quit if blank
        int cl = octant_check_list[0][icl];
        if (cl == -1) { break; }
        // Check all 8 octants to check this cluster appears in each
#pragma omp parallel for default(none) shared(dom, octant_check_list, cl) reduction(+: sum)
        for (int ioct = 0; ioct < 8; ioct++) {
            for (int jcl = 0; jcl < dom->n_clusters; jcl++) {
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

