\#!/bin/bash
#
#SBATCH --job-name=nbody_par
#SBATCH --nodes=1
#SBATCH --ntasks=1

threads=8
ts=20000
mac=3

make

echo "$ts Time Steps - Degree of $mac"
echo

echo "Serial Wall Times"
for x in 4 8 16 32 64 128; do
./nbody_psm -m $mac -Tn $ts nBodyInput$x.txt
done

echo "Serial Wall Times (GCC Optimized)":
for x in 4 8 16 32 64 128; do
./nbody_psm_opt -m $mac -Tn $ts nBodyInput$x.txt
done
