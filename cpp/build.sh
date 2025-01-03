#!/usr/bash

cat build.sh
export numThreads=10
mkdir -p build
cd build && cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH=$CONDA_PREFIX ..
make
cp src/hf ..
cp compile_commands.json ..
cd ..
./hf ./data/t1 $numThreads

