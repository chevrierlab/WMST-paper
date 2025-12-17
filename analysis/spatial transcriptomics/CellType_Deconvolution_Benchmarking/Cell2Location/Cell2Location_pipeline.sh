#!/bin/bash
#SBATCH --job-name=Cell2Loc
#SBATCH --partition=caslake
#SBATCH --account=pi-nchevrier
#SBATCH --array=0-15             
#SBATCH --cpus-per-task=11
#SBATCH --mem=165G
#SBATCH --time=36:00:00
#SBATCH --mail-user=clevenger1@rcc.uchicago.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/cell2loc_error/cell2loc_%A_%a.out
#SBATCH --error=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/cell2loc_error/cell2loc_%A_%a.err

set -euo pipefail

# Paths
ENV_DIR="/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/conda_env"
SCRIPT_DIR="/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RScripts/script"
PY_SCRIPT="${SCRIPT_DIR}/Cell2Location_pipeline.py"


# Threads hygiene
export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10
export OPENBLAS_NUM_THREADS=10

module load python/anaconda-2022.05
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
else
  echo "[ERROR] conda not found after module load"; module list; exit 1
fi

set +u
export target_platform="${target_platform:-linux-64}"
conda activate "$ENV_DIR"
echo "CONDA_PREFIX=${CONDA_PREFIX:-unset}"
export PATH=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/conda_env/bin:$PATH

unset PYTHONPATH
export PYTHONNOUSERSITE=1

# package check
python - <<'PY'
import sys, importlib
print("[DEBUG] sys.path[0] =", sys.path[0])
m = importlib.import_module("cell2location")
print("[DEBUG] cell2location file:", getattr(m, "__file__", "n/a"))
PY

echo "[SLURM] Job=${SLURM_ARRAY_JOB_ID:-NA} Task=${SLURM_ARRAY_TASK_ID:-0}"
export ARRAY_ID="${SLURM_ARRAY_TASK_ID:-0}"

# Clean stale bytecode 
find "$SCRIPT_DIR" -maxdepth 1 -name "__pycache__" -type d -exec rm -rf {} + || true
find "$SCRIPT_DIR" -maxdepth 1 -name "cell2location.*py[co]" -delete || true

# Run
python "$PY_SCRIPT"