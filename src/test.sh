
x=4
lim=2048

module load mpi

echo "Weak Scaling - 8 Processes"
while [ $x -lt $lim ]; do
echo
echo "${x} Bodies"
time mpirun -np 8 ./pnb.x < ../Input_Generation/input_files/input_${x}_psm.dat > out/pnb_4_${x}.out
let x=$((x * 2))
done

let x=4
echo
echo "Weak Scaling - 32 Processes"
while [ $x -lt $lim ]; do
echo
echo "${x} Bodies"
time mpirun -np 32 ./pnb.x < ../Input_Generation/input_files/input_${x}_psm.dat > out/pnb_32_${x}.out
let x=$((x * 2))
done

echo
echo "Strong Scaling - 512 Bodies"
for i in 1 2 4 8 16 32 64 128; do
echo
echo "${i} Processes"
time mpirun -np ${i} ./pnb.x < ../Input_Generation/input_files/input_1024_psm.dat > out/pnb_032_${i}.out
done

echo
echo "Stress Test - 128 Processes - 1024 Bodies"
time mpirun -np 128 ./pnb.x < ../Input_Generation/input_files/input_1024_psm.dat > out/pnb_stress.out
