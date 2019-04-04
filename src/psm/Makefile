default: nbody_psm nbody_psm_par

nbody_psm: nbody_psm_serial.c nbody_psm.h
	gcc --std=c99 -Wall -o nbody_psm nbody_psm_serial.c -lm

nbody_psm_par: nbody_psm_par.c nbody_psm.h
	gcc -g -O2 --std=c99 -fopenmp -Wall -o nbody_psm_parallel nbody_psm_par.c -lm

clean:
	rm -f nbody_psm
	rm -f nbody_psm_parallel
	rm -f nbody_psm_par_serial
	rm -f body*
