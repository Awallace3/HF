#!/usr/bash

echo 'Ensure you have activated a conda environment that has psi4 installed'
cat example.sh
export numThreads=10
mkdir -p build
cd build && cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH=$CONDA_PREFIX ..
make
cp src/hf ..
cp compile_commands.json ..
cd ..
cd ../psi4/
python3 ./main.py
cd ../cpp/
./hf ./data/t1 $numThreads

