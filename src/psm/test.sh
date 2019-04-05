#!/bin/bash

make

echo "Serial Timing:"
for x in 4 8 16 32 64 128; do
./nbody_psm -Tn 20000 nBodyInput$x.txt
done
