from image_processing_functions import *
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import skimage 
import os


adj = pd.read_csv("WMST_position_adjustments.csv")
adj.columns = adj.columns.str.strip()

sample_list = ["CTRL_1", "CTRL_2"]
for i in sample_list:
    
    adata = sc.read_10x_mtx(
    f'/scratch/midway2/dcipurko/STARSolo/{i}/Solo.out/Gene/raw',  
    var_names='gene_symbols',
    cache=False)
    
    tissue_mask, holes = Detect_Tissue(image_filename = f"../WM_Images/HE/{i}.png", 
                                       show_plots = True, 
                                       threshold_value = 0.84, 
                                       sigma = 2,  
                                       min_hole_size = 1500)

    
    index_pos = adj[adj['Sample'] == i].index[0]
    
    scaled_adata = Scale_coordinates(
        ST_object = adata,
        height_ST = adj.at[index_pos, "Cy5 Height"],
        width_ST = adj.at[index_pos, "Cy5 Width"],
        height_HE = adj.at[index_pos, "HE Height"],
        width_HE =  adj.at[index_pos, "HE Width"],
        HE_x = adj.at[index_pos, "HE X Pos"] ,
        HE_y = adj.at[index_pos, "HE Y Pos"] ,
        ST_x = adj.at[index_pos, "Cy5 X Pos"] ,
        ST_y = adj.at[index_pos, "Cy5 Y Pos"] ,
        path_to_barcode_positions = "Barcodes.csv"
    )
        
    scaled_filtered_adata = Filter_ST_Object(ST_object = scaled_adata, 
                     HE_filled_mask = tissue_mask, 
                     HE_holes_mask = holes,
                     path_to_region_selection_image = f"../WM_Images/Manual_Region/{i}.png")
    
    path = f"adata_objects/h5ad/WM_Raw_Scaled_Coords/"
    if not os.path.exists(path):
        os.makedirs(path)

    scaled_filtered_adata.write(f"adata_objects/h5ad/WM_Raw_Scaled_Coords/{i}.h5ad")
    
    adata_clustered = Cluster_ST_Object(ST_object = scaled_filtered_adata, resolution =1.1)
    
    path = f"adata_objects/h5ad/WM_Pre_Organ_Annotation/"
    if not os.path.exists(path):
        os.makedirs(path)    
    
    adata_clustered.write(f"adata_objects/h5ad/WM_Pre_Organ_Annotation/{i}.h5ad")
