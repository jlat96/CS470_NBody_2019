# CS 470 Research - NBody Problem
## src
This directory contains all source code for the three algorithms we are going to analyze.
### bhs
This directory contains the Barnes-Hut simulation code, which is written in C++.
### cubep3m
This directory contains the CUBEP<sup>3</sup> simulation code, which is written in Fortran 90.
### psm
This directory contains the Parker-Sochacki simulation code, which was originally written in Matlab and has been ported to C. The directory contains a makefile for compiling the serial program, and the input file 'sol_system.txt'.

NOTE: At this time, the file IO for nbody_psm_serial.c is unreliable, and hard-coded input values are used. This will be fixed before the final deliverable

Also in this folder is a copy of the code used by a group for an n-body project in 2017, which was used for reference.

### Description of Standard Program Input
An input file begins with an integer value denoting the number of bodies in the system
The next 7 values are the body's mass, the x, y, and z values for its position relative to the center of the system, and the body's u, v , and w velocity values.

A sample file with 2 bodies would appear as follows

2\
1\
0\
0\
0\
0\
0\
0\
7.6923076923076926E-9\
-12.10226300000\
-26.73256000000\
6.362842900000\
0.1714028159354\
-.1021868903979\
-.3854629379438E-01\
