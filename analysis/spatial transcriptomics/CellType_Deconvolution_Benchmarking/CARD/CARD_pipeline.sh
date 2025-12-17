#!/usr/bin/env bash
#SBATCH --job-name=CARD
#SBATCH --partition=caslake
#SBATCH --account=pi-nchevrier
#SBATCH --array=0-3             
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=36:00:00
#SBATCH --mail-user=clevenger1@rcc.uchicago.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/error/CARD_%A_%a.out
#SBATCH --error=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/error/CARD_%A_%a.err



set -euo pipefail

# Threads hygiene
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export R_FUTURE_FORK_ENABLE=false

module unload R/4.2.1
module unload openblas || true
module load openblas/0.3.13
module load udunits/2.2
module load geos/3.9.1
module load gdal/3.3.3
module load gcc/10.2.0
module load sqlite
module load R/4.2.1

GDAL_PREFIX="$(gdal-config --prefix)"
export PKG_CONFIG_PATH="${GDAL_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="${GDAL_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

gdal-config --version
pkg-config --modversion sqlite3

echo "[SLURM] Job=${SLURM_ARRAY_JOB_ID:-NA} Task=${SLURM_ARRAY_TASK_ID:-0}"

# Pass array index
export ARRAY_ID="${SLURM_ARRAY_TASK_ID:-0}"

Rscript /project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/script/CARD_pipeline.R
# Run
