#!/bin/bash

ts=10000
p=1
x=100
lim=801

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


echo "Weak Scaling"
while [ $x -lt $lim ]; do
echo
echo "${p} Processes ${x} Bodies"
time srun -n ${p} ./n_body -f ../../Input_Generation/input_files_log_10/input_${x}_bhs.in -s ${ts}
let p=$((p * 2))
let x=$((x * 2))
done


echo
echo "Stress Test - 128 Processes - 1024 Bodies"
time srun -n 128 ./n_body -f ../../Input_Generation/input_files_pow/input_1024_bhs.in -s ${ts}
