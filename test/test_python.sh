#! /bin/bash

py_exe=$(which python)
if [ -z $py_exe ]; then
    echo "Python3 not found"
    exit 1
fi
$py_exe ../python/lammpsIO.py --input ../data/dump.lammpstrj
