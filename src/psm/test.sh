\#!/bin/bash
#
#SBATCH --job-name=nbody_par
#SBATCH --nodes=1
#SBATCH --ntasks=1

threads=8
ts=50000
mac=3

make

echo
echo "$ts Time Steps - Maclaurin Degree of $mac"
echo

echo "Serial Wall Times"
for x in 4 8 16 32 64 128; do
./nbody_psm -m $mac -Tn $ts nBodyInput$x.txt
done

echo "Serial Wall Times (GCC Optimized)":
for x in 4 8 16 32 64 128; do
./nbody_psm_opt -m $mac -Tn $ts nBodyInput$x.txt
done

echo "OpenMP Serial Wall Times"
for x in 4 8 16 32 64 128; do
./nbody_psm_serial -m $mac -Tn $ts nBodyInput$x.txt
done

echo "OpenMP - Weak Scaling - $threads threads"
for x in 4 8 16 32 64 128; do
OMP_NUM_THREADS=$threads ./nbody_psm_parallel -m $mac -Tn $ts nBodyInput$x.txt
done

echo "OpenMP - Strong Scaling: "
for x in 1 2 4 8; do
OMP_NUM_THREADS=$x ./nbody_psm_parallel -m $mac -Tn $ts nBodyInput64.txt
done

