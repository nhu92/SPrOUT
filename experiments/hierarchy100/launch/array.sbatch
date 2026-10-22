#!/bin/bash
#SBATCH --job-name=h100_array
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-99
#SBATCH --output=/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022/logs/array_%A_%a.out
#SBATCH --error=/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022/logs/array_%A_%a.err
set -euo pipefail
ROOT=/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022
export PATH=/home/nhu/miniconda3/envs/sprout/bin:$PATH
export PYTHONPATH=$ROOT/code:$ROOT/sources/MrKnowItAll/src:$ROOT/sources
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python -m experiments.hierarchy100.runner run --root "$ROOT" --index "$SLURM_ARRAY_TASK_ID"
