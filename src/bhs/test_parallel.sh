cp sol_system.txt build/
cd build
mpirun -np 4 ./n_body -o positions_parallel.txt -f sol_system.txt -s 10000 -v
