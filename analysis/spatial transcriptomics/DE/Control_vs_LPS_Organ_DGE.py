import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import skimage 
import os
import gc

adata = sc.read_h5ad("adata_objects/h5ad/WM_Post_Expression_Correction/All_Samples.h5ad")

adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.pp.filter_genes(adata, min_counts = 25)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

for organ in adata.obs["Organ"].unique():
    adata_organ = adata[adata.obs["Organ"] == organ]
    
    sc.tl.rank_genes_groups( 
        adata_organ,  
        groupby="Treatment", 
        method='wilcoxon'
    ) 

    dge_df = sc.get.rank_genes_groups_df(adata_organ, group= "LPS")

    path = f"DGE_Data/Control_vs_LPS_Organ/"
    if not os.path.exists(path):
        os.makedirs(path)

    dge_df.to_csv(f"DGE_Data/Control_vs_LPS_Organ/{organ}.csv")
    del adata_organ
    gc.collect()