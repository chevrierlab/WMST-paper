# Import Libraries
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os
import gc
from multiprocessing import Pool

# Load data
adata_filtered_Celltypes = sc.read_h5ad("objects/adata_filtered_Celltypes.h5ad")

# Scanpy CPUs
sc.settings.n_jobs = -1

# Cell Type level dictionary
CellTypes = {
    "Subtype": "CellType_Subset",
    "Broad": "CellType_Broad"
}

# Format
adata_filtered_Celltypes.obs["Organ"] = adata_filtered_Celltypes.obs["Organ"].astype(str)
adata_filtered_Celltypes.obs["Subregion"] = adata_filtered_Celltypes.obs["Subregion"].astype(str)

# DE function
def process_celltype(i, Organ, Subregion, level, column, cell_type):
    print(f"Task {i + 1}: {Organ} _ {Subregion} _ {level} _ {cell_type}", flush=True)

    # Subset data
    adata_subset = adata_filtered_Celltypes[
        (adata_filtered_Celltypes.obs["Organ"] == Organ) &
        (adata_filtered_Celltypes.obs["Subregion"] == Subregion) &
        (adata_filtered_Celltypes.obs[column] == cell_type)
    ].copy()


    unique_treatments = adata_subset.obs["Treatment"].unique()
    treatment_counts = adata_subset.obs["Treatment"].value_counts()

    # Run DE  if Treatment has >1 unique value
    if len(unique_treatments) > 1 and (treatment_counts >= 2).all():

        print(f"Starting task {i + 1}: {Organ} _ {Subregion} _ {level} _ {cell_type}", flush=True)
        
        path = f"{level}/{Organ}/{Subregion}/"
        filename = f"{cell_type.replace(' ', '_')}.csv"
        full_path = os.path.join(path, filename)

        if os.path.exists(full_path):
            print(f"Previously completed task {i + 1}: {Organ} _ {Subregion} _ {level} _ {cell_type} at: {full_path}", flush=True)
            return

        adata_subset.var_names_make_unique()
        sc.pp.calculate_qc_metrics(adata_subset, qc_vars=["mt"], inplace=True)
        sc.pp.filter_genes(adata_subset, min_counts=25)

        if adata_subset.shape[1] == 0:
            print(f"No genes left after filtering for task {i + 1}: {Organ} _ {Subregion} _ {level} _ {cell_type}", flush=True)
            return  # skip

        sc.pp.normalize_total(adata_subset, inplace=True)
        sc.pp.log1p(adata_subset)

        # Perform DE analysis
        sc.tl.rank_genes_groups(adata_subset, groupby="Treatment", method="wilcoxon")
        DEG_df = sc.get.rank_genes_groups_df(adata_subset, group="LPS")

        # Define output path
        os.makedirs(path, exist_ok=True)

        # Save file with the cell type name
        DEG_df.to_csv(full_path, index=False)

        print(f"Completed task {i + 1}: {Organ} _ {Subregion} _ {level} _ {cell_type} at: {os.path.join(path, filename)}", flush=True)
    
    else:
         print(f"Did not complete task {i + 1}: {Organ} _ {Subregion} _ {level} _ {cell_type}", flush=True)
    
    # Free memory
    del adata_subset
    gc.collect()

# Prepare parallel tasks
raw_tasks = [
    (Organ, Subregion, level, column, cell_type)
    for Organ in adata_filtered_Celltypes.obs["Organ"].dropna().unique()
    for Subregion in adata_filtered_Celltypes[adata_filtered_Celltypes.obs["Organ"] == Organ].obs["Subregion"].dropna().unique()
    for level, column in CellTypes.items()
    for cell_type in adata_filtered_Celltypes[
        (adata_filtered_Celltypes.obs["Organ"] == Organ) &
        (adata_filtered_Celltypes.obs["Subregion"] == Subregion)
    ].obs[column].dropna().unique().astype(str)
]

total_tasks = len(raw_tasks)

# Add indexing to each task
tasks = [
    (i, *task) for i, task in enumerate(raw_tasks)
]

# Run parallel jobs
if __name__ == '__main__':
    with Pool(processes=47) as pool:  
        pool.starmap(process_celltype, tasks)
