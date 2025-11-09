#!/bin/bash
set -x
rm -rf build/
mkdir build && cd build
cmake -DBUILD_DOC=ON -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 ..
make -j 6
make doc_doxygen
./MolSim ../input/eingabe-collision-big.txt 5 0.0002
