#!/bin/bash

MKLROOT=/opt/intel/compilers_and_libraries/mac/mkl

LIBMKL1=$MKLROOT/lib/libmkl_intel_lp64.a
LIBMKL2=$MKLROOT/lib/libmkl_sequential.a
LIBMKL3=$MKLROOT/lib/libmkl_core.a

rm -rf macos
mkdir macos

TARGET=macos/mlm-gwas

g++ *.cpp -o $TARGET -O2 -std=c++11 $LIBMKL1 $LIBMKL2 $LIBMKL3 -lpthread -lm -ldl
