#!/usr/bin/env bash

x=2
lim=2048

while [ $x -lt $lim ]; do
    python body_generator.py --all -m 4 --output=input_files/input $x
    let x=$((x * 2))
done