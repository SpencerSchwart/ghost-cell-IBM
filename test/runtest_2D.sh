#!/bin/bash
set -e

# Compile all cases first. 

echo "Compiling sessile-ibm.c . . ."
./compile.sh sessile-ibm.c -O2 -events -display=1

echo "Compiling sessile-inclined-ibm.c . . ."
./compile.sh sessile-inclined-ibm.c -O2 -events -display=1 

echo "Compiling sessile-cylinder-ibm.c . . ."
./compile.sh sessile-cylinder-ibm.c -O2 -events -display=1

echo "Compiling sessile-sphere-ibm.c . . ."
./compile.sh sessile-sphere-ibm.c -O2 -events -display=1

# Run each case.

cd sessile-ibm
./sessile-ibm 1> out 2> log &

cd ../sessile-inclined-ibm
./sessile-inclined-ibm 1> out 2> log &

cd ../sessile-cylinder-ibm
./sessile-cylinder-ibm 1> out 2> log &

cd ../sessile-sphere-ibm
./sessile-sphere-ibm 1> out 2> log &

echo -e "\n\nWating . . .\n"
wait

echo "All Cases finished!"

