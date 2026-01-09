#!/bin/bash
# Build Script for CoolMUC (without VTK)
# Run from project root

set -x

module load slurm_setup
module load cmake/3.30.0 ninja/1.12.1 gcc/14.2.0

rm -rf build/
mkdir build && cd build || exit 1

# Configure CMake (release, no docs, no vtk)
cmake -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_DOC=OFF \
    -DENABLE_VTK_OUTPUT=OFF \
    ..

cmake --build .

echo "Build complete"
