\#!/bin/bash
#
#SBATCH --job-name=nbody_par
#SBATCH --nodes=1
#SBATCH --ntasks=1

make

echo "Serial Wall Times"
for x in 4 8 16 32 64 128; do
./nbody_psm -Tn $1 nBodyInput$x.txt
done

echo "Serial Wall Times (GCC Optimized)":
for x in 4 8 16 32 64 128; do
./nbody_psm_opt -Tn $1 nBodyInput$x.txt
done

echo "OpenMP Serial Wall Times:"
for x in 4 8 16 32 64 128; do
./nbody_psm_serial -Tn $1 nBodyInput$x.txt
done

echo "OpenMP - Weak Scaling:"
for x in 4 8 16 32 64 128; do
OMP_NUM_THREADS=8 ./nbody_psm_parallel -Tn $1 nBodyInput$x.txt
done

echo "OpenMP - Strong Scaling:"
for x in 1 2 4 8; do
OMP_NUM_THREADS=$x ./nbody_psm_parallel -Tn $1 nBodyInput64.txt
done

