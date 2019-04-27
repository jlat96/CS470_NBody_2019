#!/bin/bash

module load mpi
mpif90 -O3 -o pnb.x pnb.f90
