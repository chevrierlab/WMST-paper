import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import skimage 
import os
import gc
import matplotlib.colors as mcolors


adata = sc.read_h5ad("adata_objects/h5ad/WM_Post_Expression_Correction/All_Samples.h5ad")
adata = adata[adata.obs["Subregion"] != "Other"]
gc.collect()

for subregion in adata.obs["Subregion"].unique():
    subregion_adata = adata[adata.obs["Subregion"] == subregion]

    if len(subregion_adata.obs["Treatment"].unique()) > 1:
        subregion_adata.var_names_make_unique()
        subregion_adata.var["mt"] = subregion_adata.var_names.str.startswith("mt-")
        sc.pp.calculate_qc_metrics(subregion_adata, qc_vars=["mt"], inplace=True)
        sc.pp.filter_genes(subregion_adata, min_counts = 25)
        sc.pp.normalize_total(subregion_adata, inplace=True)
        sc.pp.log1p(subregion_adata)
        
        sc.tl.rank_genes_groups( 
            subregion_adata,  
            groupby="Treatment", 
            method='wilcoxon'
        ) 
        
        
        
        dge_df = sc.get.rank_genes_groups_df(subregion_adata, group= "LPS")

        path = f"DGE_Data/Control_vs_LPS_Subcluster_DGE/"
        if not os.path.exists(path):
            os.makedirs(path)
        if subregion == "Midbrain/Hindbrain":
            dge_df.to_csv(f"DGE_Data/Control_vs_LPS_Subcluster_DGE/Midbrain_Hindbrain.csv")     
        else:
            dge_df.to_csv(f"DGE_Data/Control_vs_LPS_Subcluster_DGE/{subregion}.csv")
            
    del subregion
    gc.collect()

