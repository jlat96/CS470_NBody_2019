#!/bin/bash

ts=10000
p=1
x=100
lim=2000

module load mpi

echo "Warning Shot"
time srun -n 8 ./n_body -f ../../Input_Generation/input_files_log_10/input_4_bhs.in -s ${ts}


echo "Barnes-Hut Timing Results"
echo
echo "Strong Scaling - 100 Bodies"
for i in 1 2 4 8 16 32 64 128; do
echo
echo "${i} Processes"
time srun -n ${i} ./n_body -f ../../Input_Generation/input_files_log_10/input_100_bhs.in -s ${ts}
done


echo "Strong Scaling - 1000 Bodies"
for i in 1 2 4 8 16 32 64 128; do
echo
echo "${i} Processes"
time srun -n ${i} ./n_body -f ../../Input_Generation/input_files_log_10/input_1000_bhs.in -s ${ts}
done
