default: nbody_psm

nbody_psm: pnb.f90
	mpif90 -O3 -o pnb.x pnb.f90

clean:
	rm pnb.x
	rm -f particle*
