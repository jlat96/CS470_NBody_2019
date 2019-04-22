#!/usr/bin/env bash

x=100
lim=1001

while [ $x -lt $lim ]; do
    python body_generator.py --all -m 4 --output=input_files_2/input $x
    let x=$x+100
done