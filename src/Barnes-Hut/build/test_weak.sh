#!/bin/bash

ts=10000

module load mpi

echo "Barnes-Hut Timing Results - 10000 Timesteps"
echo
echo "Weak Scaling - 8 Processes"
for x in 4 8 16 32 64 128 256 512 1024; do
echo
echo "${x} Bodies"
time srun -n 8 ./n_body -f ../../Input_Generation/input_files/input_${x}_bhs.in -s ${ts}
done

echo
echo "Weak Scaling - 32 Processes"
for x in 4 8 16 32 64 128 256 512 1024; do
echo
echo "${x} Bodies"
time srun -n 32 ./n_body -f ../../Input_Generation/input_files/input_${x}_bhs.in -s ${ts}
done
