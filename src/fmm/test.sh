#!/bin/bash
#
#SBATCH --job-name=par_fmm
#SBATCH --output=par_fmm_results.txt
 
echo "FMM Timing Results - 10000 Timesteps"
echo
echo "Weak Scaling - 1-8 Threads"
threads=1
for x in 100 200 400 800; do
	echo
	echo "${x} Bodies ${threads} threads"
	OMP_NUM_THREADS=${threads} salloc ./gravity N ${x} file_name ../Input_Generation/input_files_log_10/input_${x}_fmm.in  
	threads=$(expr "$threads" '*' 2)
done

echo
echo "Strong Scaling - 100 Bodies"
for x in 1 2 4 8; do
	echo
	echo "${x} Threads"
	OMP_NUM_THREADS=${x} salloc ./gravity N 100 file_name ../Input_Generation/input_files_log_10/input_100_fmm.in
done

echo
echo "Stress Test - 8 Threads - 1024 Bodies"
OMP_NUM_THREADS=8 salloc ./gravity N 1024 file_name ../Input_Generation/input_files_pow/input_1024_fmm.in


