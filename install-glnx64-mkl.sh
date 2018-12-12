#!/bin/bash

MKLROOT=/opt/intel/compilers_and_libraries/linux/mkl

LIBMKL1=$MKLROOT/lib/intel64/libmkl_intel_lp64.a
LIBMKL2=$MKLROOT/lib/intel64/libmkl_sequential.a
LIBMKL3=$MKLROOT/lib/intel64/libmkl_core.a

rm -rf glnx64
mkdir glnx64

TARGET=glnx64/mlm-gwas

g++ *.cpp -o $TARGET -s -O2 -std=c++11 -static -fopenmp \
    -Wl,--start-group $LIBMKL1 $LIBMKL2 $LIBMKL3 -Wl,--end-group \
    -lpthread -lm -ldl
