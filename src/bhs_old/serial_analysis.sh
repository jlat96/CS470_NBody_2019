#!/bin/bash

cp nBodyInput1.txt build/nBodyInput1.txt
cp nBodyInput2.txt build/nBodyInput2.txt
cp nBodyInput* build/

cd build

echo "Serial Timing:"
for x in 1 2 4 8 16 32 64 128; do
    echo "salloc -n 1 ./n_body -o outputs.txt -f nBodyInput$x.txt -s 20000 -v"
    time salloc -n 1 ./n_body -o outputs.txt -f nBodyInput$x.txt -s 20000 -v
    echo ""
    echo "================================="
done
