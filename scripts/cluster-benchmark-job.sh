#!/bin/bash
#SBATCH -J MolSim_Benchmark
#SBATCH -M cm4
#SBATCH -p cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH -o %x.%j.%N.out
#SBATCH -e %x.%j.%N.err

module load slurm_setup
module load cmake/3.30.0 ninja/1.12.1 gcc/14.2.0

cd "$SLURM_SUBMIT_DIR/build" || exit 1

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false

srun ./MolSim ../input/profiling.yaml benchmark info linked