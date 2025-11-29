#!/bin/bash
set -x
rm -rf build/
# shellcheck disable=SC2164
mkdir build && cd build
cmake -DBUILD_DOC=ON -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 ..
cmake --build .
make doc_doxygen
./MolSim ../input/eingabe-collision.txt 5 0.0002 file info P:OFF