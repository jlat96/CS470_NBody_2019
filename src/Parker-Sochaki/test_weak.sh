
x=4
lim=2048

module load mpi

echo "Parker-Sochaki Timing Results"
echo
echo "Weak Scaling - 8 Processes"
while [ $x -lt $lim ]; do
echo
echo "${x} Bodies"
time srun -n 8 ./pnb.x < ../Input_Generation/input_files/input_${x}_psm.dat > out/pnb_4_${x}.out
let x=$((x * 2))
done

let x=4
echo
echo "Weak Scaling - 32 Processes"
while [ $x -lt $lim ]; do
echo
echo "${x} Bodies"
time srun -n 32 ./pnb.x < ../Input_Generation/input_files/input_${x}_psm.dat > out/pnb_32_${x}.out
let x=$((x * 2))
done
