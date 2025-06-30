#!/bin/bash

g++ -std=c++17 -Wall -g -pedantic -I ./include  \
    src/scaleEqns.cpp                           \
    src/interpolate1d.cpp                       \
    src/hcubature.c                             \
    src/linearFit.cpp                           \
    src/bubbleScale.cpp                         \
    src/selfenergy.cpp                          \
    src/writeFile.cpp                           \
    src/curved.cpp                              \
    -o bin/Curved.elf -lgsl -lgslcblas -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -ldl

echo "Compilation finished!"
echo

#./bin/Curved.elf 399

#for (( i = 245; i>=205; i-=5));
#do
    #echo $i
    #./bin/Curved.elf $i
#done

#echo "Program executed!"



