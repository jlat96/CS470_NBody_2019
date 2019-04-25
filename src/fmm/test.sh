#!/bin/bash
#
#SBATCH --job-name=par_fmm
#SBATCH --output=par_fmm_results.txt
 
echo "FMM Timing Results - 10000 Timesteps"
echo
echo "Weak Scaling - 8 Threads"
for x in 4 8 16 32 64 128 256 512 1024; do
	echo
	echo "${x} Bodies"
	OMP_NUM_THREADS=8 salloc ./gravity N ${x} file_name ../Input_Generation/input_files_pow/input_${x}_fmm.in  
done

echo
echo "Strong Scaling - 512 Bodies"
for x in 1 2 4 8; do
	echo
	echo "${x} Threads"
	OMP_NUM_THREADS=${x} salloc ./gravity N 512 file_name ../Input_Generation/input_files_pow/input_512_fmm.in
done

echo
echo "Stress Test - 8 Threads - 1024 Bodies"
OMP_NUM_THREADS=8 salloc ./gravity N 1024 file_name ../Input_Generation/input_files_pow/input_1024_fmm.in


