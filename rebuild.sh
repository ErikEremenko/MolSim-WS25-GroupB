#!/bin/bash
set -x
rm -rf build/
mkdir build
cd build
cmake -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 ..
make -j 6
./MolSim ../input/eingabe-sonne.txt 1000 0.014
