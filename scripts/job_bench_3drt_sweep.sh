#!/bin/bash
#SBATCH -J MolSim_Bench_3DRT
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=01:00:00
#SBATCH -o ./%x.%j.%N.out
#SBATCH -e ./%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE

set -euo pipefail

# Load modules
module load slurm_setup
module load cmake/3.30.0 ninja/1.12.1 gcc/14.2.0

cd "$SLURM_SUBMIT_DIR/build" || exit 1

echo "=========================================="
echo "Submit dir: ${SLURM_SUBMIT_DIR}"
echo "Work dir:   $(pwd)"
echo "Host:       $(hostname)"
echo "Binary:     $(ls -l ./MolSim 2>/dev/null || true)"
echo "Input file: $(ls -l ../input/rayleigh_taylor_3d.yaml 2>/dev/null || true)"
echo "CPUs/task:  ${SLURM_CPUS_PER_TASK}"
echo "=========================================="

if [[ ! -x ./MolSim ]]; then
  echo "ERROR: ./MolSim not found or not executable in $(pwd)" >&2
  exit 2
fi

if [[ ! -f ../input/rayleigh_taylor_3d.yaml ]]; then
  echo "ERROR: ../input/rayleigh_taylor_3d.yaml not found (relative to $(pwd))" >&2
  exit 3
fi

# Pinning / reproducibility
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DYNAMIC=false

# Required thread counts
THREADS_LIST=(1 2 4 8 14 16 28 56)
STRATEGIES=("S:COLORING" "S:TASKBASED")

OUTDIR="output/benchmarks"
mkdir -p "$OUTDIR"
CSV="$OUTDIR/bench_3drt_${SLURM_JOB_ID}.csv"

echo "job_id,host,strategy,omp_threads,elapsed_s,total_iterations,mups" > "$CSV"

run_one() {
  local strategy="$1"
  local threads="$2"

  export OMP_NUM_THREADS="$threads"

  local tag="${strategy}_${threads}T"
  local log="$OUTDIR/bench_3drt_${SLURM_JOB_ID}_${tag}.log"

  echo "------------------------------------------"
  echo "Strategy: ${strategy}"
  echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
  echo "Log: ${log}"
  echo "------------------------------------------"

  # YAML mode usage:
  # ./MolSim filename [file|benchmark] [loglevel] [linked|direct] [P:ON|P:OFF] [S:COLORING|S:TASKBASED]
  # explicitly override parallelization+strategy to avoid YAML drift.
  srun --cpu-bind=cores ./MolSim ../input/rayleigh_taylor_3d.yaml benchmark info linked P:ON "$strategy" 2>&1 | tee "$log"

  local elapsed iters mups
  elapsed=$(grep -m1 "Time elapsed:" "$log" | awk '{print $(NF-1)}')
  iters=$(grep -m1 "Total iterations:" "$log" | awk '{print $NF}')
  mups=$(grep -m1 "Molecule-Updates per Second" "$log" | awk '{print $NF}')

  echo "${SLURM_JOB_ID},$(hostname),${strategy},${threads},${elapsed},${iters},${mups}" >> "$CSV"
}

for strategy in "${STRATEGIES[@]}"; do
  for threads in "${THREADS_LIST[@]}"; do
    run_one "$strategy" "$threads"
  done
done

echo "=========================================="
echo "Benchmark sweep complete"
echo "CSV: ${CSV}"
echo "Logs: ${OUTDIR}/bench_3drt_${SLURM_JOB_ID}_*.log"
echo "=========================================="
