module load mpi

lim=801
p=1
x=100

echo "Strong Scaling - 100 Bodies"
for i in 1 2 4 8 16; do
    echo
    echo "${i} Processes"
    time mpirun -np ${i} ./build/n_body -f ../Input_Generation/input_files_log_10/input_100_bhs.in -s 10000
done

echo
echo "Weak Scaling"
while [ $x -lt $lim ]; do
    echo
    echo "${p} Processes ${x} Bodies"
    time mpirun -np ${p} ./build/n_body -f ../Input_Generation/input_files_log_10/input_${x}_bhs.in -s 10000
let p=$((p * 2))
let x=$((x * 2))
done
