#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "types.h"
#include "comms.h"
#include "input.h"
#include "output.h"

void comms_initialise(int* argc, char*** argv, MpiInfo* info) {
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &info->n_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &info->my_rank);
    DT_ARGS = create_Args_datatype();
    DT_DECOMPRESULTS = create_DecompResults_datatype();
    DT_RUNRESULTS = create_RunResults_datatype();
}

Args comms_read_input(int argc, char** argv, MpiInfo* info) {
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

    if (info->my_rank == 0) {
        // Read command line arguments
        int ra = read_args(argc, argv, &args);
        if (ra == 0) {
            printf("# Command line arguments successfully processed\n");
        } else {
            printf("Error Processing command line arguments (error code %i)\n", ra);
            abort();
        }
    }

    if (MPI_Bcast(&args, 1, DT_ARGS, 0, MPI_COMM_WORLD)) {
        fprintf(stderr, "Failure broadcasting arguments.\n");
        abort();
    }

    return args;
}

// no longer using a cartesian communicator; instead each process handles its own run
//void comms_cell_map(CellDomain* dom, MpiInfo* info) {
//    int ndims = 3;
//    int reorder = 0;
//    // no periodic boundary conditions
//    int pbc[3] = {0, 0, 0};
//
//    int proot = (int) ceil(cbrt((double) info->n_processors));
//    int dims[3] = {proot, proot, proot};
//
//    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, pbc, reorder, &info->cart_comm);
//    for (int i = -1; i <= 1; i++) {
//        for (int j = -1; j <= 1; j++) {
//            for (int k = -1; k <= 1; k++) {
//                int rank = info->my_rank;
//                MPI_Cart_shift(info->cart_comm, 0, i, &rank, &rank);
//                MPI_Cart_shift(info->cart_comm, 1, j, &rank, &rank);
//                MPI_Cart_shift(info->cart_comm, 2, k, &rank, &rank);
//                info->neighbours[i][j][k] = rank;
//            }
//        }
//    }
//
//    info->n_routers = dom->cell_n_routers[info->my_rank];
//    info->routers = dom->cell_routers[info->my_rank];
//}

void comms_finalise() {
    MPI_Type_free(&DT_RUNRESULTS);
    MPI_Type_free(&DT_DECOMPRESULTS);
    MPI_Type_free(&DT_ARGS);
    MPI_Finalize();
}

MPI_Datatype create_Args_datatype() {
    int blocklengths[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint displacements[8] = {
            offsetof(Args, router_radius),
            offsetof(Args, space_station_initial_size),
            offsetof(Args, space_station_size_increment),
            offsetof(Args, number_space_station_size_increments),
            offsetof(Args, volume_fraction_initial),
            offsetof(Args, volume_fraction_increment),
            offsetof(Args, number_volume_fraction_increments),
            offsetof(Args, use_current_time_as_seed),
    };
    MPI_Datatype datatypes[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT,
                                 MPI_C_BOOL};
    MPI_Datatype type;
    MPI_Type_create_struct(8, blocklengths, displacements, datatypes, &type);

    MPI_Aint lb, extent;
    MPI_Type_get_extent(type, &lb, &extent);
    MPI_Type_create_resized(type, lb, extent, &type);
    MPI_Type_commit(&type);
    return type;
}

MPI_Datatype create_DecompResults_datatype() {
    int blocklengths[6] = {1, 1, 1, 1, 1, 1};
    MPI_Aint displacements[6] = {
            offsetof(DecompResults, ncells),
            offsetof(DecompResults, cluster_search_time),
            offsetof(DecompResults, cluster_merge_time),
            offsetof(DecompResults, cluster_span_time),
            offsetof(DecompResults, n_clusters),
            offsetof(DecompResults, spanning_cluster)
    };
    MPI_Datatype datatypes[6] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
    MPI_Datatype type;
    MPI_Type_create_struct(6, blocklengths, displacements, datatypes, &type);

    MPI_Aint lb, extent;
    MPI_Type_get_extent(type, &lb, &extent);
    MPI_Type_create_resized(type, lb, extent, &type);
    MPI_Type_commit(&type);
    return type;
}

MPI_Datatype create_RunResults_datatype() {
    // Send the value of the pointer over as bytes, but it will be overwritten by the receiving code
    int blocklengths[4] = {1, 1, 1, sizeof(DecompResults*)};
    MPI_Aint displacements[4] = {
            offsetof(RunResults, S),
            offsetof(RunResults, P),
            offsetof(RunResults, n_cell_sizes),
            offsetof(RunResults, decomp_results)
    };
    MPI_Datatype datatypes[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_BYTE};
    MPI_Datatype type;
    MPI_Type_create_struct(4, blocklengths, displacements, datatypes, &type);

    MPI_Aint lb, extent;
    MPI_Type_get_extent(type, &lb, &extent);
    MPI_Type_create_resized(type, lb, extent, &type);
    MPI_Type_commit(&type);
    return type;
}
