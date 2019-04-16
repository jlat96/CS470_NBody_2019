This is a parallel implementation of the Parker-Sochacki algorithm in C ported from Matlab. This directory also contains an example of the format we will be using for our input data.

## Building/Installation
Use the included Makefile to build

## Running
Usage: ./nbody_psm_serial [-dghntTv]\
	-d: Print debug output\
	-g: Specify the granularity of the output (how often to report new state)\
	-h: Print usage\
	-n: Specify the number of time steps to complete\
	-t: Specify the time step\
	-T: Print timer output\
	-v: Print verbose output

