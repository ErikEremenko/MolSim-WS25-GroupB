#!/bin/bash
set -x
rm -rf build/
# shellcheck disable=SC2164
mkdir build && cd build
cmake -DBUILD_DOC=OFF -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 ..
make -j 6
./MolSim ../input/eingabe-collision.txt 5 0.0002 benchmark off