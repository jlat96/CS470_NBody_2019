#!/bin/bash
#
# Project/Account (change to your own)
#SBATCH -A hpc2n-1234-56
#
# Number of nodes
#SBATCH -N 8
#
# Use nodes exclusive
#SBATCH --exclusive
#
#

# Load the compiler and MPI library you compiled the program with. Here, openmpi/gcc   
module load mpi

# Total number of MPI tasks will be calculated by slurm based on either the defaults or command line parameters.

srun ./n_body -f ../../Input_Generation/input_files/input_512_bhs.in -s ${ts}


# End of submit file
