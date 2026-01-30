#!/bin/bash
set -x
rm -rf build/
# shellcheck disable=SC2164
mkdir build && cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_DOC=ON -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 ..
cmake --build .
cmake --build . --target doc_doxygen
# sets CPUs to maximum available frequency
sudo cpupower frequency-set -g performance
# taskset -c binds the process to the specified set of CPU cores
# Using single-threaded calculation
sudo taskset -c 19 chrt -r 50 nice -n -10 \
./MolSim ../input/argon_equilibrate_smoothed.yaml file info linked
# reset the ownership of all files in build to the current user
cd .. && sudo chown -R "$USER":"$USER" build/