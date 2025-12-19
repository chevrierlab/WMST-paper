#!/bin/bash
#SBATCH --job-name=pyscenic
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=16
#SBATCH --mem=90G
#SBATCH --time=36:00:00
#SBATCH --account=pi-nchevrier
#SBATCH --mail-user=clevenger1@rcc.uchicago.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/PythonScripts/error/pyscenic_%A_%a.out
#SBATCH --error=/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/PythonScripts/error/pyscenic_%A_%a.err


# CONFIG
ENV="/project/nchevrier/projects/clevenger1/software/conda_envs/pyscenic"
WHITELIST="/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/whitelist/tfs.txt"
DB10="/project/nchevrier/projects/clevenger1/software/conda_envs/pyscenic_reference/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
DB05="/project/nchevrier/projects/clevenger1/software/conda_envs/pyscenic_reference/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
ANN="/project/nchevrier/projects/clevenger1/software/conda_envs/pyscenic_reference/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
IN_DIR="/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/Celltypes/Objects/core"
OUT="/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/Celltypes/output/core"
NWORKERS="${NWORKERS:-16}"

mkdir -p "$OUT/logs"

# Thread hygiene
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Numba avoid TBB
export NUMBA_THREADING_LAYER=workqueue
export NUMBA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Enb
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate "$ENV"
echo "[info] using pyscenic: $(which pyscenic)"

shopt -s nullglob
for LPS in "$IN_DIR"/*_LPS.h5ad; do
  base="$(basename "$LPS" .h5ad)"     
  ORG="${base%_LPS}"                  
  ODIR="$OUT/$ORG"
  mkdir -p "$ODIR"

  if [[ ! -s "$LPS" ]]; then
    echo "[skip] Missing LPS file for $ORG: $LPS"
    continue
  fi

  echo "[run] ${ORG} | File: $LPS"
  # GRN
  pyscenic grn "$LPS" "$WHITELIST" \
    --method grnboost2 \
    --output "$ODIR/${ORG}.adjacencies.tsv" \
    --num_workers "$NWORKERS" \
    2>&1 | tee "$OUT/logs/${ORG}_grn.log"

  # CTX 
  pyscenic ctx "$ODIR/${ORG}.adjacencies.tsv" \
    "$DB10" "$DB05" \
    --annotations "$ANN" \
    --expression_mtx "$LPS" \
    --output "$ODIR/${ORG}.regulons.csv" \
    --num_workers "$NWORKERS" \
    2>&1 | tee "$OUT/logs/${ORG}_ctx.log"


done


