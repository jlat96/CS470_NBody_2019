#!/bin/bash

make

echo "Serial Timing:"
for x in 1 2 4 8 16 32 64 128; do
./nbody_psm -T nBodyInput$1.txt
done
