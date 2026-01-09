#!/bin/bash
set -x
# sets CPUs to maximum available frequency
sudo cpupower frequency-set -g performance
# taskset -c binds the process to the specified set of CPU cores
# Using single-threaded calculation
sudo taskset -c 0 chrt -r 50 nice -n -10 \
perf record -g -o perf.data -- ./MolSim ../input/profiling.yaml benchmark info linked