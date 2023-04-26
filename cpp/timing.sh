#!/usr/bash

# rm -rf build
# rm hf
# mkdir -p build
# cd build && cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug ..
# make
# cp src/hf ..
# cp compile_commands.json ..
# cd ..

export dp="data/"
export dn="t1"
export fn=$dp$dn

echo "\nTiming for $fn\n" > outputs/$dn.out
for i in 1 2 4 6 8 10
do
    export OMP_NUM_THREADS=$i
    echo ./hf $fn $i >> outputs/$dn.out
    ./hf $fn $i >> outputs/$dn.out
done

export dn="t3"
export fn=$dp$dn
echo "Timing for $fn" > outputs/$dn.out
for i in 1 2 4 6 8 10
do
    export OMP_NUM_THREADS=$i
    echo ./hf $fn $i >> outputs/$dn.out
    ./hf $fn $i >> outputs/$dn.out
done
