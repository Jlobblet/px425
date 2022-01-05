#ifndef PX425_TYPES_H
#define PX425_TYPES_H

#include <stdbool.h>
#include "mpi.h"

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
    /// Size and length along each dimnesion
    double S, Lx, Ly, Lz;
    /// Number of cells along each dimension and sizes of cells
    int nx, ny, nz, nc;
    /// Routers in each cell
    Router*** cell_routers;
    /// Number of routers in each cell
    int* cell_n_routers;
    /// Total number of routers
    int n_routers;
    /// Cluster number
    int n_clusters;
    /// True if clusters span the space
    bool spanning_cluster;
} CellDomain;

typedef struct MpiInfo {
    int n_processors;
    int my_rank;
} MpiInfo;

void generate_routers(int* Nrtr, Router** Rtr, double S, double R, double P);

void create_domain_decomp(CellDomain* dom, int Nrtr, Router* Rtr, int nx, int ny, int nz);

void destroy_domain_decomp(CellDomain* dom);

#endif //PX425_TYPES_H
