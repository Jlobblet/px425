#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4571mb
#SBATCH --time=00:02:00
#SBATCH --reservation PX425

# NOTE: For MPI, we change ntasks-per-node and keep  cpus-per-task at 1

# Change into the directory where I submitted the job
cd "$SLURM_SUBMIT_DIR" || exit 1

# Load the modules used when compiling my code
module purge
module load intel/2017.2.174-GCC-5.3.0-2.27 impi/2017.2.174

# Note we use srun rather than mpirun here
srun ./wifi -S 150 -N 1 -P 0.345838 -M 3 -Q 0.000001

exit 0
