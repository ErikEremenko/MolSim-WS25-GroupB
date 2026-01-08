#!/bin/bash
set -x
rm -rf build/
# shellcheck disable=SC2164
mkdir build && cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_DOC=ON -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 ..
cmake --build .