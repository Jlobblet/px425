//  Comms routines for PX425 assignment 4. Contains all
//  routines which interact with MPI libraries. Many of
//  these are currently incomplete and will work only in
//  serial. You will need to correct this.
//
//  Original code created by N. Hine  - November 2021
//  (based on previous code by D. Quigley)
#include "mpi.h"
#include "grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "comms.h"

// Number of processors
int p;
// Rank of current processor
int my_rank;

// Cartesian communicator
MPI_Comm cart_comm;

// Coordinates of current rank in the processor grid
int my_rank_coords[2];

// Ranks of neighbours to the current processor (left, right, down, up)
int my_rank_neighbours[4];

// Time and initialisation and shutdown
double t1, t2;

/// Function to initialise MPI, get the communicator size p and the
/// rank my_rank of the current process within that communicator.
/// N. Hine (based on code by D. Quigley) - Univ. Warwick
void comms_initialise(int argc, char** argv) {
    // square root of p
    int proot;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // Start the timer. Set t1 using MPI_Wtime() which returns a double
    t1 = MPI_Wtime();

    // Check that we have a square number of processors
    proot = (int) sqrt((double) p + 0.5);
    if (proot * proot != p) {
        if (my_rank == 0) {
            printf("Number of processors must be an exact square!\n");
            exit(EXIT_FAILURE);
        }
    }
}

/// Function to map our p processors into a 2D Cartesian grid of
/// dimension proot by proot where proot = sqrt(p).
///
/// Should populate the arrays my_rank_cooords, which contains the
/// location of the current MPI task within the processor grid, and
/// my_rank_neighbours, which contains (in the order left, right,
/// down and up) ranks of neighbouring MPI tasks on the grid with
/// which the current task will need to communicate.
/// N. Hine (based on code by D. Quigley) - University of Warwick
void comms_processor_map() {
    // Information for setting up a Cartesian communicator
    int ndims = 2;
    int reorder = 0;
    int pbc[2] = {1, 1};
    int dims[2];
    int my_cart_rank;

    // Local variables
    // square root of p
    int proot;
    // Square root of number of processors
    proot = (int) sqrt((double) p + 0.5);

    // Dimensions of Cartesian communicator
    dims[x] = proot;
    dims[y] = proot;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, pbc, reorder, &cart_comm);
    MPI_Comm_rank(cart_comm, &my_cart_rank);
    MPI_Cart_coords(cart_comm, my_cart_rank, 2, my_rank_coords);

    MPI_Cart_shift(cart_comm, 0, -1, my_rank_coords, &my_rank_neighbours[left]);
    MPI_Cart_shift(cart_comm, 0, 1, my_rank_coords, &my_rank_neighbours[right]);
    MPI_Cart_shift(cart_comm, 1, -1, my_rank_coords, &my_rank_neighbours[down]);
    MPI_Cart_shift(cart_comm, 1, 1, my_rank_coords, &my_rank_neighbours[up]);

//    printf("My rank       : %d\n"
//           "My coordinates: %d %d\n"
//           "Up            : %d\n"
//           "Down          : %d\n"
//           "Left          : %d\n"
//           "Right         : %d\n"
//           "\n",
//           my_cart_rank,
//           my_rank_coords[0],
//           my_rank_coords[1],
//           my_rank_neighbours[up],
//           my_rank_neighbours[down],
//           my_rank_neighbours[left],
//           my_rank_neighbours[right]
//           );
}

/// Function to compute the glocal magnetisation of the grid by
/// averaging over all values of local_mag, and storing the result
/// in global_mag.
/// N. Hine (based on code by D. Quigley) - University of Warwick
void comms_get_global_mag(double local_mag, double* global_mag) {
    // This is only correct on one processor. You will need
    // to use a collective communication routine to correct this.
//    *global_mag = local_mag;


    // Insert collective communication operation here
    MPI_Allreduce(&local_mag, global_mag, 1, MPI_DOUBLE, MPI_SUM, cart_comm);


    *global_mag = *global_mag / (double) p;

}

/// Function to send boundary spins on each side of the local grid
/// to neighbour processors, and to receive from those processors the
/// halo information needed to perform computations involving spins
/// on the boundary processors grid.
/// N. Hine (based on code by D. Quigley) - University of Warwick
void comms_halo_swaps() {
    // Send and receive buffers
    int* sendbuf, * recvbuf;

    // MPI Status
    MPI_Status status;

    // Loop counters
    int ix, iy;

    // If running on 1 processor copy boundary elements into opposite halo
    if (p == 1) {

        for (iy = 0; iy < grid_domain_size; iy++) {
            grid_halo[right][iy] = grid_spin[iy][0];
            grid_halo[left][iy] = grid_spin[iy][grid_domain_size - 1];
        }

        for (ix = 0; ix < grid_domain_size; ix++) {
            grid_halo[up][ix] = grid_spin[0][ix];
            grid_halo[down][ix] = grid_spin[grid_domain_size - 1][ix];
        }

        // Do not do any comms
        return;
    }

    // Allocate buffers
    sendbuf = (int*) malloc(grid_domain_size * sizeof(int));
    if (sendbuf == NULL) {
        printf("Error allocating sendbuf in comms_halo_swaps\n");
        exit(EXIT_FAILURE);
    }
    recvbuf = (int*) malloc(grid_domain_size * sizeof(int));
    if (recvbuf == NULL) {
        printf("Error allocating recvbuf in comms_halo_swaps\n");
        exit(EXIT_FAILURE);
    }

    // SASHAY LEFT
    for (iy = 0; iy < grid_domain_size; iy++) {
        // Copy the left-hand elements into sendbuf
        sendbuf[iy] = grid_spin[iy][0];
    }
    MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], left,
                 recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], left,
                 cart_comm, &status);
    // Copy recvbuf into the right-hand elements
    memcpy(grid_halo[right], recvbuf, grid_domain_size * sizeof(int));

    // SHIMMY RIGHT
    for (iy = 0; iy < grid_domain_size; iy++) {
        // Copy the right-hand elements into sendbuf
        sendbuf[iy] = grid_spin[iy][grid_domain_size - 1];
    }
    MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], right,
                 recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], right,
                 cart_comm, &status);
    memcpy(grid_halo[left], recvbuf, grid_domain_size * sizeof(int));

    // PRANCE FORWARD
    // Since the grid is stored this way around, can memcpy it.
    memcpy(sendbuf, grid_spin[0], grid_domain_size * sizeof(int));
    MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], down,
                 recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], down,
                 cart_comm, &status);
    memcpy(grid_halo[up], recvbuf, grid_domain_size * sizeof(int));

    // BOOGIE DOWN
    memcpy(sendbuf, grid_spin[grid_domain_size - 1], grid_domain_size * sizeof(int));
    MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], up,
                 recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], up,
                 cart_comm, &status);
    memcpy(grid_halo[down], recvbuf, grid_domain_size * sizeof(int));

    // Release memory
    free(sendbuf);
    free(recvbuf);
}


/// Function to collect all contributions to the global grid onto
/// rank zero for visualisation.
/// N. Hine (based on code by D. Quigley) - University of Warwick
void comms_get_global_grid() {
    // comms buffer
    int* combuff;

    // MPI Status
    MPI_Status status;

    // Information on the remote domain
    int remote_domain_start[2] = {0, 0};

    // Loop counters and error flag
    int ix, iy, ixg, iyg, ip;

    // Just point at local grid if running on one processor
    if (p == 1) {
        global_grid_spin = grid_spin;
        return;
    }

    if (my_rank == 0) {

        // Rank 0 first fills out its part of the global grid
        for (iy = 0; iy < grid_domain_size; iy++) {
            for (ix = 0; ix < grid_domain_size; ix++) {

                // Global indices
                ixg = ix + grid_domain_start[x];
                iyg = iy + grid_domain_start[y];

                global_grid_spin[iyg][ixg] = grid_spin[iy][ix];

            }
        }
    }

    // Allocate buffer
    combuff = (int*) malloc(grid_domain_size * sizeof(int));
    if (combuff == NULL) {
        printf("Error allocating combuff in comms_get_global_grid\n");
        exit(EXIT_FAILURE);
    }

    if (my_rank == 0) {

        // Now loops over all other ranks receiving their data
        for (ip = 1; ip < p; ip++) {

            // First receive remote_domain_start from rank ip
            MPI_Recv(remote_domain_start, 2, MPI_INT, ip, MPI_ANY_TAG, cart_comm, &status);

            // Loop over rows within a domain
            for (iy = 0; iy < grid_domain_size; iy++) {

                // Receive this row from rank ip
                MPI_Recv(combuff, grid_domain_size, MPI_INT, ip, MPI_ANY_TAG, cart_comm, &status);

                for (ix = 0; ix < grid_domain_size; ix++) {
                    // Global indices
                    ixg = ix + remote_domain_start[x];
                    iyg = iy + remote_domain_start[y];

                    // Store in global_grid_spin
                    global_grid_spin[iyg][ixg] = combuff[ix];
                } // elements in row
            } // rows
        } // processors

    } else {
        // All other processors must send the data rank 0 needs
        // Send grid_domain_start to rank 0
        MPI_Send(grid_domain_start, 2, MPI_INT, 0, 0, cart_comm);

        // Loop over rows in the domain, sending them to rank 0
        for (iy = 0; iy < grid_domain_size; iy++) {
            memcpy(combuff, grid_spin[iy], grid_domain_size * sizeof(int));
            MPI_Send(combuff, grid_domain_size, MPI_INT, 0, 0, cart_comm);
        }
    }

    free(combuff);
}


/// Function to finalise MPI functionality and exit cleanly
/// N. Hine (based on code by D. Quigley) - University of Warwick
void comms_finalise() {
    // Measure the time t2 using MPI_Wtime() which returns a double
    t2 = MPI_Wtime();

    if (my_rank == 0 && p > 1) {
        printf("Total time elapsed since MPI initialised :  %12.6f s\n", t2 - t1);
    }

    MPI_Finalize();
}
