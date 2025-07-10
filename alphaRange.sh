#!/bin/bash

g++ -std=c++17 -Wall -g -pedantic -I ./include  \
    src/scaleEqns.cpp                           \
    src/interpolate1d.cpp                       \
    src/hcubature.c                             \
    src/linearFit.cpp                           \
    src/bubbleScale.cpp                         \
    src/selfenergy.cpp                          \
    src/radialScale.cpp                         \
    src/seRad.cpp                               \
    src/tangentialScale.cpp                     \
    src/seTang.cpp                              \
    src/writeFile.cpp                           \
    src/alphaRange.cpp                           \
    -o bin/alphaRange.elf -lgsl -lgslcblas -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -ldl

echo "Compilation finished!"
echo

./bin/alphaRange.elf

echo "Program executed!"



