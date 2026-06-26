#!/bin/bash
set -e

# Compile all cases first. All cases are in serial except for the inclined plates

echo "Compiling sessile-ibm3D.c . . ."
./compile.sh sessile-ibm3D.c -O2 -events -display=1

echo "Compiling sessile-inclined-ibm3D.c . . ."
./compile.sh sessile-inclined-ibm3D.c -O2 -events -display=1 -mpi

echo "Compiling sessile-sphere-ibm3D.c . . ."
./compile.sh sessile-sphere-ibm3D.c -O2 -events -display=1

# Run each case.

cd sessile-ibm3D
./sessile-ibm3D 1> out 2> log &

cd ../sessile-inclined-ibm3D
mpirun -n 3 ./sessile-inclined-ibm3D 1> out 2> log &

cd ../sessile-sphere-ibm3D
./sessile-sphere-ibm3D 1> out 2> log &

echo "Wating . . ."
wait

echo "All Cases finished!"

