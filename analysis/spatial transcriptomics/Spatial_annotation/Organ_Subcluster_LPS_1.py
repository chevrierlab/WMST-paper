import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import skimage 
import os
import gc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


for sample in ["LPS_1"]:
    adata = sc.read_h5ad(f"adata_objects/h5ad/WM_Post_Expression_Correction/{sample}.h5ad")
    adata = adata[adata.obs["Organ"] != "Other"]

    image = skimage.io.imread(f"../WM_Images/HE/{sample}.png")[:,:,:3]
    img = skimage.color.rgb2gray(image)

    for organ in adata.obs["Organ"].unique():
        adata_sub = adata[adata.obs["Organ"] == organ]

        if organ == "SK":
            adata_sub = adata_sub[adata_sub.obs["Selected"] == True]
            
        adata_sub.var_names_make_unique()
        adata_sub.var["mt"] = adata_sub.var_names.str.startswith("mt-")
        sc.pp.calculate_qc_metrics(adata_sub, qc_vars=["mt"], inplace=True)
        sc.pp.filter_genes(adata_sub, min_counts = 25)
        sc.pp.normalize_total(adata_sub, inplace=True)
        sc.pp.log1p(adata_sub)
        sc.pp.highly_variable_genes(adata_sub, flavor="seurat", n_top_genes=3000)

        sc.pp.pca(adata_sub)
        sc.pp.neighbors(adata_sub)
        print("Running UMAP")
        sc.tl.umap(adata_sub)
        if organ in ["LI", "SK"]:
            sc.tl.leiden(adata_sub, resolution = 0.8, key_added="clusters")
        else:
            sc.tl.leiden(adata_sub, resolution = 1, key_added="clusters")
        sc.pl.umap(adata_sub, color="clusters", show=False)
        path = f"plots/clustering/Single_Organ/{sample}/"
        if not os.path.exists(path):
            os.makedirs(path)
        plt.savefig(f"plots/clustering/Single_Organ/{sample}/UMAP_{organ}_Clusters.png", bbox_inches='tight')

        adata_sub_meta = adata_sub.obs

        x_bottom_percentile, x_top_percentile = np.percentile(adata_sub_meta["x_scaled_image"], [4, 96])
        y_bottom_percentile, y_top_percentile = np.percentile(adata_sub_meta["y_scaled_image"], [4, 96])

        x_range = (x_top_percentile - x_bottom_percentile)*1.8
        y_range = (y_top_percentile - y_bottom_percentile)*1.8

        max_dim = int(np.max([x_range, y_range])/2)
        x_middle = int(np.mean([x_bottom_percentile, x_top_percentile]))
        y_middle = int(np.mean([y_bottom_percentile, y_top_percentile]))

        y_lower = max(0, y_middle-max_dim)
        y_upper = min(len(img), y_middle+max_dim)

        x_lower = max(0, x_middle-max_dim)
        x_upper = min(len(image[1]), x_middle+max_dim)


        img_cropped = img[y_lower:y_upper, x_lower:x_upper]
        HE_image_cropped = image[y_lower:y_upper, x_lower:x_upper]

        if max_dim > 3500: 
            img_cropped = img
            HE_image_cropped = image
            adata_sub_meta["x_scaled_image"] =adata_sub_meta["x_scaled_image"]
            adata_sub_meta["y_scaled_image"] = adata_sub_meta["y_scaled_image"]
            low, high = np.percentile(img_cropped, (0.1, 60))
            img_cropped = skimage.exposure.rescale_intensity(img_cropped, in_range=(low, high))
        else:
            adata_sub_meta["x_scaled_image"] = adata_sub_meta["x_scaled_image"] - x_lower
            adata_sub_meta["y_scaled_image"] = adata_sub_meta["y_scaled_image"] - y_lower
            low, high = np.percentile(img_cropped, (0.1, 95))
            img_cropped = skimage.exposure.rescale_intensity(img_cropped, in_range=(low, high))


        min_cropped_dim = min(len(img_cropped), len(img_cropped[1]))
        point_size = 1.93**((len(image[1])) / min_cropped_dim) - 0.6

        if  min_cropped_dim < 350:
            point_size = 60
        if min_cropped_dim >= 350 and min_cropped_dim <= 425:
            point_size = 40
        if min_cropped_dim >= 425 and min_cropped_dim < 650:
            point_size = 22


        sns.relplot(
            data = adata_sub_meta, 
            x = "x_scaled_image", 
            y = "y_scaled_image", 
            hue = "clusters", 
            linewidth=0,
            s = point_size, 
            height=19, 
            aspect=0.4
        )


        plt.imshow(img_cropped, cmap="gray")
        plt.axis("off")

        plt.savefig(f"plots/clustering/Single_Organ/{sample}/{organ}_Clusters.png", bbox_inches='tight')

        path = f"../WM_Images/Single_Organ/{sample}/HE/"
        if not os.path.exists(path):
            os.makedirs(path)
        skimage.io.imsave(fname=f"../WM_Images/Single_Organ/{sample}/HE/{organ}.png", arr=HE_image_cropped)

        path = f"../WM_Images/Single_Organ/{sample}/BW/"
        if not os.path.exists(path):
            os.makedirs(path)
        skimage.io.imsave(fname=f"../WM_Images/Single_Organ/{sample}/BW/{organ}.png", arr=img_cropped)


        adata_sub.obs = adata_sub_meta
        path = f"adata_objects/h5ad/Single_Organ/{sample}/"
        if not os.path.exists(path):
            os.makedirs(path)
        adata_sub.write(f"adata_objects/h5ad/Single_Organ/{sample}/{organ}.h5ad")


        #DEG analysis:

        sc.tl.rank_genes_groups( 
            adata_sub,  
            groupby="clusters", 
            method='wilcoxon'
        ) 

        for k in adata_sub.obs["clusters"].unique():
            dge_df = sc.get.rank_genes_groups_df(adata_sub, group= k)

            path = f"DGE_Data/Organ_subclusters/{sample}/{organ}/"
            if not os.path.exists(path):
                os.makedirs(path)

            dge_df.to_csv(f'DGE_Data/Organ_subclusters/{sample}/{organ}/Cluster_{k}_DGE.csv')

            for h in dge_df["names"].tolist()[0:15]: 
                adata_gene_index = adata_sub.var.index.get_loc(h)
                m = adata_sub.X[:,adata_gene_index].max()
                adata_sub.obs["gene"] = np.concatenate(adata_sub.X[:,adata_gene_index].toarray().tolist()).flat
                alpha = adata_sub.X[:,adata_gene_index].toarray().tolist()/m

                path = f"plots/DGE/Organ_Subclusters/{sample}/{organ}/Cluster_{k}"
                if not os.path.exists(path):
                    os.makedirs(path)

                sns.relplot(
                    data = adata_sub.obs, 
                    x = "x_scaled_image", 
                    y = "y_scaled_image", 
                    hue = "gene", 
                    linewidth=0,
                    s = point_size, 
                    alpha = alpha,
                    palette = sns.color_palette("magma", as_cmap=True), 
                    legend = False
                )

                plt.imshow(img_cropped, cmap="gray")
                plt.axis("off")


                plt.savefig(f"plots/DGE/Organ_Subclusters/{sample}/{organ}/Cluster_{k}/{h}.png", bbox_inches='tight')

                plt.close()
                plt.clf()
                gc.collect()

        del adata_sub
        gc.collect()

