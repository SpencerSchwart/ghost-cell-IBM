#!/bin/bash
set -e

# Compile all cases first. 

#echo "Compiling impact1 . . ."
#./compile.sh impact-sphere.c impact1 -O2 -events -display=1 -mpi -g

echo "Compiling impact2 . . ."
./compile.sh impact-sphere.c impact2 -O2 -events -display=1 -mpi -g

echo "Compiling impact3 . . ."
./compile.sh impact-sphere.c impact3 -O2 -events -display=1 -mpi -g

# Run each case.

# cd impact1
# ./run.sh -mpi 3 &
# tail -f out
# echo -e "\n\nImpact 1 done.\n"

cd impact2
./run.sh -mpi 3
echo -e "\n\nImpact 2 done.\n"

cd ../impact3
./run.sh -mpi 4
echo -e "\n\nImpact 3 done.\n"

echo "All Cases finished!"

