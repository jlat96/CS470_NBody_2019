module load mpi

np=$1

time srun -n ${np} mpirun ./pnb.x < pnb_032.dat > pnb_032.out &
