module load mpi

np=$1

time mpirun -np $np ./pnb.x < pnb_032.dat > pnb_032.out &
