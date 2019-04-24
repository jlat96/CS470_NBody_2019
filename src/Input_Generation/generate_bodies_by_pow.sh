#!/usr/bin/env bash

x=2
lim=2049

while [ $x -lt $lim ]; do
    python body_generator.py --all -m 4 --output=input_files_pow/input $x
    let x=$((x * 2))
done
