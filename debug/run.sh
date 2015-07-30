#!/bin/bash

# After using this script it is necessary to run again the buil.sh script
# for generating again the *.mod files

source ../config.sh

DEBUGFLAGS="-g -fbounds-check -O0"

\rm -f exe
\rm -f *.o
\rm -f ../src/*.o
\rm -f ../src/*.mod

$FF $DEBUGFLAGS -c ../src/*.f90
$FF $DEBUGFLAGS -I$(pwd)/../src/ -c *.f90
$FF $DEBUGFLAGS -I../src/ -o exe *.o

\rm -f ../src/*.o
\rm -f ../src/*.mod

./exe
