#!/bin/bash
#SBATCH -J MolSim_Bench_3DRT_1
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=00:30:00
#SBATCH -o ./%x.%j.%N.out
#SBATCH -e ./%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE

set -euo pipefail

module load slurm_setup
module load cmake/3.30.0 ninja/1.12.1 gcc/14.2.0

cd "$SLURM_SUBMIT_DIR/build" || exit 1

if [[ ! -x ./MolSim ]]; then
  echo "ERROR: ./MolSim not found or not executable in $(pwd)" >&2
  exit 2
fi

if [[ ! -f ../input/rayleigh_taylor_3d.yaml ]]; then
  echo "ERROR: ../input/rayleigh_taylor_3d.yaml not found (relative to $(pwd))" >&2
  exit 3
fi

# Choose thread count + strategy (override via sbatch --export)
# Example:
#   sbatch --cpus-per-task=14 --export=ALL,THREADS=14,STRATEGY=S:COLORING scripts/job_bench_3drt_single.sh
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK}}
STRATEGY=${STRATEGY:-S:COLORING}

export OMP_NUM_THREADS="$THREADS"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DYNAMIC=false

OUTDIR="output/benchmarks"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/bench_3drt_${SLURM_JOB_ID}_${STRATEGY}_${OMP_NUM_THREADS}T.log"

echo "=========================================="
echo "Submit dir: ${SLURM_SUBMIT_DIR}"
echo "Work dir:   $(pwd)"
echo "Host:       $(hostname)"
echo "CPUs/task:  ${SLURM_CPUS_PER_TASK}"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
echo "Strategy:   ${STRATEGY}"
echo "Log:        ${LOG}"
echo "=========================================="

srun --cpu-bind=cores ./MolSim ../input/rayleigh_taylor_3d.yaml benchmark info linked P:ON "$STRATEGY" 2>&1 | tee "$LOG"

elapsed=$(grep -m1 "Time elapsed:" "$LOG" | awk '{print $(NF-1)}')
iters=$(grep -m1 "Total iterations:" "$LOG" | awk '{print $NF}')
mups=$(grep -m1 "Molecule-Updates per Second" "$LOG" | awk '{print $NF}')
echo "RESULT threads=${OMP_NUM_THREADS} strategy=${STRATEGY} elapsed_s=${elapsed} iters=${iters} mups=${mups}"
