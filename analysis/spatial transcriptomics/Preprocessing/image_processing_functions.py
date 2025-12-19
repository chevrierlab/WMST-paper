import scipy as sp
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import skimage
from math import floor

def Detect_Tissue(image_filename, 
                sigma = 1, 
                threshold_value = "auto", 
                connectivity = 2, 
                min_size = "auto",
                min_hole_size = False,
                show_plots = True):
    
    image = skimage.io.imread(image_filename)[:,:,:3]
    gray_image = skimage.color.rgb2gray(image)
    blurred_image = skimage.filters.gaussian(gray_image, sigma=sigma)
    
    if threshold_value == "auto":
        t = skimage.filters.threshold_otsu(blurred_image)
    elif type(threshold_value) == float:
        t = threshold_value
        
    binary_mask = blurred_image < t
    
    if min_size == "auto":
        min_size = image.size/10
    
    object_mask = skimage.morphology.remove_small_objects(binary_mask,min_size)
    holes_filled = sp.ndimage.morphology.binary_fill_holes(object_mask)
    
    if min_hole_size != False:
        holes = holes_filled.astype(np.int8) - object_mask.astype(np.int8)
        
        if min_hole_size == "auto":
            s = holes_filled.size/1500
            
        if  type(min_hole_size) == int:
            s = min_hole_size
        
        selected_holes = skimage.morphology.remove_small_objects(holes > 0, min_size = s)

        tissue_holes_removed = holes_filled.astype(np.int8) - selected_holes.astype(np.int8)
        
        if show_plots == True:
            fig = plt.figure()

            fig.add_subplot(1,3,1)
            plt.imshow(image)
            plt.title("H&E", fontsize=10)
            plt.axis('off')

            fig.add_subplot(1,3,2)
            plt.imshow(holes_filled)
            plt.title("Detected Tissue \n(Filled in)", fontsize=10)
            plt.axis('off')           

            fig.add_subplot(1,3,3)
            plt.imshow(tissue_holes_removed)
            plt.title("Detected Tissue \n(Holes Removed)", fontsize=10)
            plt.axis('off')
            
            plt.tight_layout()
            
        return holes_filled, selected_holes
    
    else:
        if show_plots == True:
        
            plt.subplot(1,2,1)
            plt.imshow(image)
            plt.title("H&E")
            plt.axis('off')

            plt.subplot(1,2,2)
            plt.imshow(holes_filled)
            plt.title("Detected Tissue (Filled in)")
            plt.axis('off')
            plt.show()
            
            plt.tight_layout()
            
        return holes_filled



def Scale_coordinates(ST_object, 
                    height_HE, 
                    width_HE, 
                    height_ST,
                    width_ST, 
                    HE_x, 
                    HE_y, 
                    ST_x, 
                    ST_y, 
                    path_to_barcode_positions):
    
    
    HE_position_x = round(HE_x - (width_HE/2))
    HE_position_y = round(HE_y - (height_HE/2))
    ST_position_x = round(ST_x - (width_ST/2))
    ST_position_y = round(ST_y - (height_ST/2))
    

    ST_shift_x = (ST_position_x - HE_position_x)
    ST_shift_y = (ST_position_y - HE_position_y)
    coords = pd.read_csv(path_to_barcode_positions)
    coords["x_scaled_image"] = coords["x_image"].multiply(width_ST/coords["x_image"].max()).add(ST_shift_x)
    coords["y_scaled_image"] = coords["y_image"].multiply(height_ST/coords["y_image"].max()).add(ST_shift_y)
    ST_object.obs = ST_object.obs.join(coords.set_index("barcodes"))
    
    return ST_object


def under_tissue_function(x, y, HE_Width, HE_Height, HE_filled_mask, blue_cannel = None, path_to_region_selection_image = None):
    if x >=0 and y >= 0 and x < HE_Width and y < HE_Height: 
        y_pixel = int(floor(y))
        x_pixel = int(floor(x))
        truth_value = HE_filled_mask[y_pixel, x_pixel] == True
        if type(path_to_region_selection_image) == str:
            if blue_cannel[y_pixel, x_pixel] > 150:
                truth_value = False
    else:
        truth_value = False
            
    return truth_value



def cavity_function(x, y, HE_Width, HE_Height, red_channel, path_to_region_selection_image, HE_holes_mask):
    if x >=0 and y >= 0 and x <= HE_Width and y <= HE_Height: 
        y_pixel = int(round(y))
        x_pixel = int(round(x))
        truth_value = HE_holes_mask[y_pixel, x_pixel] == True
        if type(path_to_region_selection_image) == str:
            if red_channel[y_pixel, x_pixel] > 150:
                truth_value = True
    else:
        truth_value = False
            
    return truth_value


def region_selection_function(x, y, HE_Width, HE_Height, green_channel, path_to_region_selection_image):
    if type(path_to_region_selection_image) == str:
        if x >=0 and y >= 0 and x <= HE_Width and y <= HE_Height: 
            y_pixel = int(round(y))
            x_pixel = int(round(x))
            if green_channel[y_pixel, x_pixel] > 150:
                truth_value = True
            else:
                truth_value = False
    else:
        truth_value = False
            
    return truth_value

def Filter_ST_Object(ST_object, 
                    HE_filled_mask, 
                    HE_holes_mask,
                    path_to_region_selection_image = None,
                    return_raw=False):
    
    HE_Height = len(HE_filled_mask)
    HE_Width = len(HE_filled_mask[1])
    raw_obj = ST_object.copy()
    cols_to_drop = [col for col in raw_obj.obs.columns if "Cavity" in col]
    raw_obj.obs = raw_obj.obs.drop(cols_to_drop, axis=1)
    raw_metadata = raw_obj.obs.copy()

    if type(path_to_region_selection_image) == str:
        selection_image = skimage.io.imread(path_to_region_selection_image)[:,:,:3]
        red_channel = selection_image[:, :, 0]
        green_channel = selection_image[:, :, 1]
        blue_channel = selection_image[:, :, 2]

        cavity_bool_list = [cavity_function(x, y, HE_Width, HE_Height, red_channel, path_to_region_selection_image, HE_holes_mask) for x, y in zip(ST_object.obs['x_scaled_image'], ST_object.obs['y_scaled_image'])]
        ST_object.obs["Cavity"] = cavity_bool_list
        print("Number of spots in a cavity region = " + str(ST_object.obs["Cavity"].sum()))

        selected_bool_list = [region_selection_function(x, y, HE_Width, HE_Height, green_channel, path_to_region_selection_image) for x, y in zip(ST_object.obs['x_scaled_image'], ST_object.obs['y_scaled_image'])]
        ST_object.obs["Selected"] = selected_bool_list
        
        under_tissue_bool_list = [under_tissue_function(x, y, HE_Width, HE_Height, blue_channel = blue_channel, path_to_region_selection_image = path_to_region_selection_image, HE_filled_mask = HE_filled_mask) for x, y in zip(ST_object.obs['x_scaled_image'], ST_object.obs['y_scaled_image'])]
        ST_object = ST_object[under_tissue_bool_list]
        print("Number of spots under tissue = " + str(len(ST_object.obs)))
        
    else: 
        under_tissue_bool_list = [under_tissue_function(x, y, HE_Width, HE_Height, HE_filled_mask = HE_filled_mask) for x, y in zip(ST_object.obs['x_scaled_image'], ST_object.obs['y_scaled_image'])]
        ST_object = ST_object[under_tissue_bool_list]
        print("Number of spots under tissue = " + str(len(ST_object.obs)))

        if "Cavity" not in ST_object.obs.columns:
            ST_object.obs["Cavity"] = False 
    
    # Perform diffusion calculations
    
    masked_metadata = ST_object.obs.copy()  # Metadata from the masked object
    print(f"Columns in masked_metadata: {list(masked_metadata.columns)}")


    # Merge Cavity column into the raw object
    merged_metadata = raw_metadata.merge(
        masked_metadata[["Barcode", "Cavity"]],
        on="Barcode",
        how="left"
    )
    print(f"Columns in merged_metadata: {list(merged_metadata.columns)}")

    merged_metadata["Cavity"] = merged_metadata["Cavity"].fillna(True)
    raw_obj.obs = merged_metadata

    # Calculate average UMI/spot in vs. out of tissue
    average_umi_per_spot = raw_obj.obs.loc[raw_obj.obs["UMI"] > 0, "UMI"].mean()
    average_umi_in_tissue_noDropouts = raw_obj.obs.loc[
        (raw_obj.obs["Cavity"] == False) & (raw_obj.obs["UMI"] > 0), 
        "UMI"
    ].mean()
    average_umi_in_tissue = raw_obj.obs.loc[
        (raw_obj.obs["Cavity"] == False), 
        "UMI"
    ].mean()
    average_umi_out_of_tissue = raw_obj.obs.loc[
        (raw_obj.obs["Cavity"] == True) & (raw_obj.obs["UMI"] > 0), 
        "UMI"
    ].mean()

    # Calculate total diffusion percent
    total_transcripts_in_cavity = raw_obj.obs.query("Cavity == True")["UMI"].sum()
    total_transcripts_under_tissue = raw_obj.obs["UMI"].sum()
    diffusion_percent = (total_transcripts_in_cavity / total_transcripts_under_tissue) * 100

    # Print formatted diffusion results with rounding
    print(f"Average UMI/spot: {average_umi_per_spot:.2f}")
    print(f"Average UMI/spot in tissue: {average_umi_in_tissue:.2f}")
    print(f"Average UMI/spot in tissue not including empty spots: {average_umi_in_tissue_noDropouts:.2f}")
    print(f"Average UMI/spot out of tissue: {average_umi_out_of_tissue:.2f}")
    print(f"Total transcripts in cavity: {total_transcripts_in_cavity:.2f}")
    print(f"Total transcripts under tissue: {total_transcripts_under_tissue:.2f}")
    print(f"Diffusion Percentage (by transcripts): {diffusion_percent:.2f}%")


    if return_raw:
        return raw_obj, ST_object  # Return both raw and filtered objects
    else:
        return ST_object  # Return only the filtered object


def Cluster_ST_Object(ST_object, 
                    min_spot_umi_count = 120, 
                    min_gene_counts = 20, 
                    n_top_genes_clustering = 7500):
    
    ST_object.var_names_make_unique()
    ST_object.var["mt"] = ST_object.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(ST_object, qc_vars=["mt"], inplace=True)
    sc.pp.filter_genes(ST_object, min_counts = min_gene_counts)

    ST_object = ST_object[(ST_object.obs["Cavity"] == False) & (ST_object.obs["total_counts"] > min_spot_umi_count)]

    sc.pp.normalize_total(ST_object, inplace=True)
    sc.pp.log1p(ST_object)
    sc.pp.highly_variable_genes(ST_object, flavor="seurat", n_top_genes=n_top_genes_clustering)

    sc.pp.pca(ST_object)
    sc.pp.neighbors(ST_object)
    print("Running UMAP")
    sc.tl.umap(ST_object)
    sc.tl.leiden(ST_object, key_added="clusters", resolution = 1.1)
    
    return ST_object
