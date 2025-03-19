#! /bin/bash

pushd ../python
    uv sync
    source .venv/bin/activate
popd
python ../python/lammpsio.py --input ../data/dump.lammpstrj
