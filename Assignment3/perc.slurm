#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2679
#SBATCH -p px425
#SBATCH --time=5:00:00
#SBATCH --job-name=perc

# The above requests a group of 4 cores, which will be
# all on one node and will share memory.

# Set OMP_NUM_THREADS. The total number of threads
# launched by your program must not exceed the
# number of cores requested with ppn above
# (be careful when using nested parallelism).
export OMP_NUM_THREADS=4

# Run the executable
./perc

# Indicate normal completion
exit 0

