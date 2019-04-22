
p=1
x=100
lim=801

module load mpi

echo "Warning Shot"
time srun -n 8 ./pnb.x < ../Input_Generation/input_files/input_4_psm.dat > out/pnb_test.out


echo "Parker-Sochaki Timing Results"
echo

echo
echo "Strong Scaling - 400 Bodies"
for i in 2 4 8 16 32 64 128; do
echo
echo "${i} Processes"
time srun -n ${i} ./pnb.x < ../Input_Generation/input_files_2/input_100_psm.dat > out/pnb_100_${i}.out 
done


echo "Weak Scaling"
while [ $x -lt $lim ]; do
echo
echo "${p} Processes ${x} Bodies"
time srun -n ${p} ./pnb.x < ../Input_Generation/input_files_2/input_${x}_psm.dat > out/pnb_4_${x}.out
let p=$((p * 2))
let x=$((x * 2))
done

echo
echo "Stress Test - 128 Processes - 1024 Bodies"
time srun -n 128 ./pnb.x < ../Input_Generation/input_files_2/input_1024_psm.dat > out/pnb_stress.out
