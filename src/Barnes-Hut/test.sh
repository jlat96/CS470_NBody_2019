
ts=10000

module load mpi

echo "Barnes-Hut Timing Results - 10000 Timesteps"
echo
echo "Weak Scaling - 8 Processes"
for x in 4 8 16 32 64 128 256 512 1024; do
echo
echo "${x} Bodies"
time mpirun -np 8 ./n_body -f ../../Input_Generation/input_files/input_${x}_bhs.in -s ${ts}
done

echo
echo "Weak Scaling - 32 Processes"
for x in 4 8 16 32 64 128 256 512 1024; do
echo
echo "${x} Bodies"
time mpirun -np 32 ./n_body -f ../../Input_Generation/input_files/input_${x}_bhs.in -s ${ts}
done

echo
echo "Strong Scaling - 512 Bodies"
for x in 4 8 16 32 64 128; do
echo
echo "${x} Processes"
time mpirun -np ${x} ./n_body -f ../../Input_Generation/input_files/input_512_bhs.in -s ${ts}
done

echo
echo "Stress Test - "Stress Test - 128 Processes - 1024 Bodies"
time mpirun -np 128 ./n_body -f ../../Input_Generation/input_files/input_1024_bhs.in -s ${ts}

