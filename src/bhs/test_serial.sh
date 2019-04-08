cp sol_system.txt build/
cd build
mpirun -np 1 ./n_body -o positions_serial.txt -f sol_system.txt -s 10000 -v
