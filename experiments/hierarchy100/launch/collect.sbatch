#!/bin/bash
#SBATCH --job-name=h100_collect
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022/logs/collect_%A_%a.out
#SBATCH --error=/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022/logs/collect_%A_%a.err
set -euo pipefail
ROOT=/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022
export PATH=/home/nhu/miniconda3/envs/sprout/bin:$PATH
export PYTHONPATH=$ROOT/code:$ROOT/sources/MrKnowItAll/src:$ROOT/sources
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python -m experiments.hierarchy100.runner collect --root "$ROOT"
