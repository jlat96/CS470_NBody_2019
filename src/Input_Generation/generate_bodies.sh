#!/usr/bin/env bash

x=2
lim=2048

while [ $x -lt $lim ]; do
    python body_generator.py --all --output=input_files/input $x
    let x=$((x * 2))
done