
p=1
x=100
lim=801

module load mpi

echo "Warning Shot"
time srun -n 8 ./pnb.x < ../Input_Generation/input_files_pow/input_4_psm.dat > out/pnb_test.out


echo "Parker-Sochaki Timing Results"
echo

echo
echo "Strong Scaling - 100 Bodies"
for i in 2 4 8 16 32 64 128; do
echo
echo "${i} Processes"
time srun -n ${i} ./pnb.x < ../Input_Generation/input_files_log_10/input_100_psm.dat > out/pnb_100_${i}.out 
done

