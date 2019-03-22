cp sol_system.txt build/
cd build
mpirun -np 1 ./n_body -o positions.txt -f sol_system.txt
