#!/bin/bash
set -x

rm -rf build/
# shellcheck disable=SC2164
mkdir build && cd build

cmake -G Ninja DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-g" -DBUILD_DOC=ON -DENABLE_VTK_OUTPUT=OFF -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 ..
cmake --build .
# sets CPUs to maximum available frequency
sudo cpupower frequency-set -g performance
# taskset -c binds the process to the specified set of CPU cores
# Using single-threaded calculation
sudo taskset -c 0 chrt -r 50 nice -n -10 \
 perf record -g --call-graph dwarf  ./MolSim ../input/rayleigh_taylor_big.yaml benchmark info linked
# reset the ownership of all files in build to the current user
cd .. && sudo chown -R "$USER":"$USER" build/