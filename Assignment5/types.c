#include "mt19937ar.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "types.h"

/// Set up number of routers for a given volume fraction, and
/// generate random coordinates for each one, setting cluster
/// value to 0 to indicate cluster not found yet.
void generate_routers(int* Nrtr, Router** Rtr, double S, double R, double P) {
    *Nrtr = (int) (S * S * S * P / (R * R * R));
    *Rtr = calloc(*Nrtr, sizeof(struct Router));
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
    dom->cell_n_routers = calloc(dom->nc, sizeof(int));
    dom->cell_routers = calloc(dom->nc, sizeof(Router*));
    dom->n_routers = Nrtr;
    dom->spanning_cluster = 0;

    // count the routers associated with each cell
    for (int i = 0; i < Nrtr; i++) {
        int ix = floor(dom->nx * Rtr[i].x / dom->Lx);
        int iy = floor(dom->ny * Rtr[i].y / dom->Ly);
        int iz = floor(dom->nz * Rtr[i].z / dom->Lz);
        int ic = ix + dom->ny * iy + dom->ny * dom->nz * iz;
        if (ic >= dom->nc) {
            printf("oh no\n");
        }
        dom->cell_n_routers[ic]++;
    }

    // allocate memory for storage of each cell's set of pointers to routers,
    // then re-set the counts of numbers of routers to zero
    for (int ic = 0; ic < dom->nc; ic++) {
        dom->cell_routers[ic] = calloc(dom->cell_n_routers[ic], sizeof(Router*));
        dom->cell_n_routers[ic] = 0;
    }

    // fill up the list of routers associated with each cell, counting them
    // again as we go for indexing purposes
    for (int i = 0; i < Nrtr; i++) {
        int ix = floor(dom->nx * Rtr[i].x / dom->Lx);
        int iy = floor(dom->ny * Rtr[i].y / dom->Ly);
        int iz = floor(dom->nz * Rtr[i].z / dom->Lz);
        int ic = ix + dom->ny * iy + dom->ny * dom->nz * iz;
        dom->cell_routers[ic][dom->cell_n_routers[ic]] = &Rtr[i];
        dom->cell_n_routers[ic]++;
    }
}

/// Deallocate storage associated with a domain decomposition struct.
void destroy_domain_decomp(CellDomain* dom) {
    for (int ic = 0; ic < dom->nc; ic++) {
        free(dom->cell_routers[ic]);
    }
    free(dom->cell_routers);
    free(dom->cell_n_routers);
}
