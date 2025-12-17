#!/bin/bash
#SBATCH --job-name=RCTD_array
#SBATCH --partition=bigmem
#SBATCH --account=pi-nchevrier
#SBATCH --array=0-3             
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=30:00:00
#SBATCH --mail-user=clevenger1@rcc.uchicago.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/error/RCTD_array_%A_%a.out
#SBATCH --error=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/error/RCTD_array_%A_%a.err


set -euo pipefail

# Threads hygiene
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export R_FUTURE_FORK_ENABLE=false

module load R/4.2.1

echo "[SLURM] Job=${SLURM_ARRAY_JOB_ID:-NA} Task=${SLURM_ARRAY_TASK_ID:-0}"

# Pass array index
export ARRAY_ID="${SLURM_ARRAY_TASK_ID:-0}"


# Run
Rscript /project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/script/RCTD_pipeline.R