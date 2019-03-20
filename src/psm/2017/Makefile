default: nbody_psm_serial

nbody_psm_serial: NBodyPSMe.c CauchyProd.h CauchyProdPower.h CauchyProdSqrng.h Horner.h NBPSMparameters.h
	gcc --std=c99 -Wall  -o nbody_psm_serial NBodyPSMe.c -lm

clean:
	rm -f nbody_psm_serial