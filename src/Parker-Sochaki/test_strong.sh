module load mpi

echo "Parker-Sochaki Timing Results"
echo

echo
echo "Strong Scaling - 512 Bodies"
for i in 1 2 4 8 16 32 64 128; do
echo
echo "${i} Processes"
time srun -n ${i} ./pnb.x < ../Input_Generation/input_files/input_1024_psm.dat > out/pnb_032_${i}.out
done
