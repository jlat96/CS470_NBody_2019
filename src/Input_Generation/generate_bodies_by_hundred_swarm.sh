#!/usr/bin/env bash

x=100
lim=2000

while [ $x -lt $lim ]; do
    python body_generator.py --all -s -m 4 --output=input_files_swarm/input $x
    let x=$x+100
done
