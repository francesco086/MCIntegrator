#!/bin/bash

source config.sh
LIBNAME="mcintegrator"

\rm -f *.a
cd src
   \rm -f *.o *.a *.mod
   $FF $OPTFLAGS -c *.f90
   ar rcv lib${LIBNAME}.a *.o
   ranlib lib${LIBNAME}.a
   mv lib${LIBNAME}.a ../
cd ..

echo
echo "Library ready!"
echo
echo "Help, how can I use it? Two options:"
echo "1)   $FF -I$(pwd)/src/ -c example.f90"
echo "     $FF -L$(pwd) example.o -l${LIBNAME}" 
echo "2)   $FF -I$(pwd)/src/ -L$(pwd) example.f90 -l${LIBNAME}"
