#!/bin/bash
#
#SBATCH --job-name=par_fmm
#SBATCH --output=par_fmm_results.txt
 
for i in 1 2 4 8; do
	echo "== test par_fmm with $i threads=="
	OMP_NUM_THREADS="$i" salloc ./gravity
done
