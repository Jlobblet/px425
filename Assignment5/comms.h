#ifndef PX425_COMMS_H
#define PX425_COMMS_H

#include "input.h"
#include "types.h"

MPI_Datatype MPI_ARGS, MPI_DECOMPRESULTS, MPI_RUNRESULTS;

MPI_Datatype create_Args_datatype();

MPI_Datatype create_DecompResults_datatype();

MPI_Datatype create_RunResults_datatype(int n_cell_sizes);

MPI_Datatype create_Results_datatype(int n_runs, MPI_Datatype MPI_RUNRESULTS);

void comms_initialise(int* argc, char*** argv, MpiInfo* info);

Args comms_read_input(int argc, char** argv, MpiInfo* info);

void comms_cell_map();

void comms_finalise();

#endif //PX425_COMMS_H
