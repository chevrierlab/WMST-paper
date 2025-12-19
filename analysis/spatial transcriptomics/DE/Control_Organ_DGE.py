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

adata = adata[adata.obs["Treatment"] == "Control"]
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.pp.filter_genes(adata, min_counts = 25)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

sc.tl.rank_genes_groups( 
    adata,  
    groupby="Organ", 
    method='wilcoxon', 
    pts = True
) 


WM_images = {}
for img_tmp in adata.obs["Sample"].unique():
    image = skimage.io.imread(f"../WM_Images/HE/{img_tmp}.png")[:,:,:3]
    image = skimage.color.rgb2gray(image)
    low, high = np.percentile(image, (0.5, 70))
    image = skimage.exposure.rescale_intensity(image, in_range=(low, high))
    WM_images[img_tmp] = image    

for organ in adata.obs["Organ"].unique():

    dge_df = sc.get.rank_genes_groups_df(adata, group= organ)

    path = f"DGE_Data/Control_Only/"
    if not os.path.exists(path):
        os.makedirs(path)

    dge_df.to_csv(f"DGE_Data/Control_Only/{organ}.csv")
    
    dge_df_organ_sub_genes = dge_df["names"][0:20]
    
    for gene in dge_df_organ_sub_genes:
        adata_gene_index = adata.var.index.get_loc(gene)
        m = adata.X[:,adata_gene_index].max()
        adata.obs["gene"] = np.concatenate(adata.X[:,adata_gene_index].toarray().tolist()).flat
        alpha = adata.X[:,adata_gene_index].toarray().tolist()/m

        for sample in adata.obs["Sample"].unique():
            adata_sample_sub = adata[adata.obs["Sample"] == sample]
            alpha_sample = alpha[adata.obs["Sample"] == sample]


            sns.relplot(
                data = adata_sample_sub.obs, 
                x = "x_scaled_image", 
                y = "y_scaled_image", 
                hue = "gene",
                linewidth=0,
                s = 1.5, 
                color='red', 
                palette=['red'],
                alpha = alpha_sample,
                legend = False,  
                height = 19, 
                aspect=0.4
            )

            plt.imshow(WM_images[sample], cmap="gray")
            plt.axis("off")

            path = f"plots/DGE/Control_Organ_Markers/{organ}/"
            if not os.path.exists(path):
                os.makedirs(path)
            plt.savefig(f"plots/DGE/Control_Organ_Markers/{organ}/{gene}_{sample}.png", bbox_inches='tight')
            plt.close()
            plt.clf()

            del adata_sample_sub
            gc.collect()
    