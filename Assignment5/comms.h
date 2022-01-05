#ifndef PX425_COMMS_H
#define PX425_COMMS_H

#include "input.h"
#include "types.h"

MPI_Datatype DT_ARGS, DT_DECOMPRESULTS, DT_RUNRESULTS;

MPI_Datatype create_Args_datatype();

MPI_Datatype create_DecompResults_datatype();

MPI_Datatype create_RunResults_datatype();

void comms_initialise(int* argc, char*** argv, MpiInfo* info);

Args comms_read_input(int argc, char** argv, MpiInfo* info);

void comms_cell_map();

void comms_finalise();

#endif //PX425_COMMS_H
