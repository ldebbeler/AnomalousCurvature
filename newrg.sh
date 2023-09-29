#!/bin/bash

g++ -std=c++17 -Wall -g -pedantic -I ./include  \
    src/scaleEqns.cpp                           \
    src/interpolate1d.cpp                       \
    src/hcubature.c                             \
    src/linearFit.cpp                           \
    src/bubbleScale.cpp                         \
    src/newSe.cpp                               \
    src/radialScale.cpp                         \
    src/seRad.cpp                               \
    src/tangentialScale.cpp                     \
    src/seTang.cpp                              \
    src/writeNew.cpp                            \
    src/newrg.cpp                               \
    -o bin/Newrg.elf -lgsl -lgslcblas -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -ldl

echo "Compilation finished!"
echo

#./bin/Newrg.elf 400
#for (( i = 400; i>=205; i-=10));
#do
#    echo $i
#    ./bin/Newrg.elf $i
#done

echo "Program executed!"



