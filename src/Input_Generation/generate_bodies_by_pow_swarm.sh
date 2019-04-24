#!/usr/bin/env bash

x=2
lim=2049

while [ $x -lt $lim ]; do
    python body_generator.py --all -s -m 4 --output=input_files_swarm/input $x
    let x=$((x * 2))
done
