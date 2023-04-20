#!/usr/bash

rm -rf build
mkdir -p build
cd build && cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug ..
make
cp src/hf ..
cp compile_commands.json ..
cd ..
./hf
