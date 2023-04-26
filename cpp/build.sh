#!/usr/bash

export numThreads=10
rm -rf build
rm hf
mkdir -p build
cd build && cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug ..
make
cp src/hf ..
cp compile_commands.json ..
cd ..
# export OMP_NUM_THREADS=$numThreads
# ./hf data/t1 $numThreads
