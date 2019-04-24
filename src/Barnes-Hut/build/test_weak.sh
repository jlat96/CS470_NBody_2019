#!/bin/bash

ts=10000
p=1
x=100
lim=2000

module load mpi

echo "Warning Shot"
time srun -n 8 ./n_body -f ../../Input_Generation/input_files_log_10/input_4_bhs.in -s ${ts}

echo "Weak Scaling"
while [ $x -lt $lim ]; do
echo
echo "${p} Processes ${x} Bodies"
time srun -n ${p} ./n_body -f ../../Input_Generation/input_files_log_10/input_${x}_bhs.in -s ${ts}
let p=$((p * 2))
let x=$((x * 2))
done

