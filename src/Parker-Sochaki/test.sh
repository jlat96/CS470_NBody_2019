#!/bin/bash

p=1
x=100
lim=801

module load mpi

echo "Warning Shot"
time srun -n 2 ./pnb.x < ../Input_Generation/input_files_pow/input_4_psm.dat > out/pnb_test.
rm comp_cost_*


echo "Parker-Sochaki Timing Results"
echo
echo "Strong Scaling - 100 Bodies"
for i in 1 2 4 8 16 32 64 128; do
echo
echo "${i} Processes"
time srun -n ${i} ./pnb.x < ../Input_Generation/input_files_log_10/input_100_psm.dat > out/pnb_100_${i}.out
rm comp_cost_*
done


echo "Weak Scaling"
while [ $x -lt $lim ]; do
echo
echo "${p} Processes ${x} Bodies"
time srun -n ${p} ./pnb.x < ../Input_Generation/input_files_log_10/input_${x}_psm.dat > out/pnb_4_${x}.out
rm comp_cost_*
let p=$((p * 2))
let x=$((x * 2))
done
