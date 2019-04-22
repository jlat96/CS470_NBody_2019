#!/bin/bash

ts=10000

module load mpi

echo "Barnes-Hut Timing Results - 10000 Timesteps"

echo
echo "Strong Scaling - 512 Bodies"
for x in 1 2 4 8 16 32 64 128; do
echo
echo "${x} Processes"
time srun -n ${x} ./n_body -f ../../Input_Generation/input_files/input_512_bhs.in -s ${ts}
done
