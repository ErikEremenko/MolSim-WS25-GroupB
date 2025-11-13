#!/bin/bash
set -x
sudo -v

rm -rf build/
# shellcheck disable=SC2164
mkdir build && cd build

cmake -DBUILD_DOC=OFF -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 ..
make -j 6
# sets CPUs to maximum available frequency
sudo cpupower frequency-set -g performance
# taskset -c binds the process to the specified set of CPU cores
# set -c argument in accordance to your machine's available core count
# Using single-threaded calculation
sudo taskset -c 0-19 chrt -r 50 nice -n -10 \
  ./MolSim ../input/eingabe-collision.txt 5 0.0002 benchmark off P:OFF