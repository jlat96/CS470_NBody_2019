
ts=20000

echo "Strong Scaling - Swarm"
for x in 4 8 16 32 64 128 256 512 1024; do
echo
echo "${x} Bodies"
time mpirun -np 32 ./n_body -f ../../Input_Generation/input_files/input_swarm_${x}_bhs.in -s ${ts}
done

echo
echo "Weak Scaling - Swarm"
for x in 4 8 16 32; do
echo
echo "${x} Threads"
time mpirun -np ${x} ./n_body -f ../../Input_Generation/input_files/input_swarm_512_bhs.in -s ${ts}
done
