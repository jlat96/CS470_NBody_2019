
p=1
x=100
lim=801

module load mpi

echo "Warning Shot"
time srun -n 8 ./pnb.x < ../Input_Generation/input_files_pow/input_4_psm.dat > out/pnb_test.out


echo "Parker-Sochaki Timing Results"
echo

echo "Weak Scaling"
while [ $x -lt $lim ]; do
echo
echo "${p} Processes ${x} Bodies"
time srun -n ${p} ./pnb.x < ../Input_Generation/input_files_log_10/input_${x}_psm.dat > out/pnb_4_${x}.out
let p=$((p * 2))
let x=$((x * 2))
done
