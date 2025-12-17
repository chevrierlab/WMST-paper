# SETUP -----------------------------------------------------

source("/project/nchevrier/projects/clevenger1/Projects/Resources/Functions_Sourcing/R_Master_Source.R")

organs = c("SP", "LU", "LI", "KI")
organ  = organs[as.integer(Sys.getenv("ARRAY_ID")) + 1]

print(organ)


# RCTD Pipeline ------------------

Run_RCTD  = function(
    organ,
    organ_dir = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Organ_exports" %>% spath,
    base_dir  = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Reference_data" %>% spath,
    rctd_dir = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RCTD_output" %>% spath,
    sample
){
  
  out_dir_check = file.path(rctd_dir, organ, paste0(organ, "_", sample, "_RCTD_proportions.csv"))
  
  if(file.exists(out_dir_check)){
    print(paste0(organ, "|", sample, ": already ran, output exists, skipping."))
    return(NULL)
  }
  
  # Load and process components
  
  comp_dir = file.path(base_dir, organ, "components")
  counts_ref = Matrix::readMM(gzfile(file.path(comp_dir, paste0(organ, "_X_counts.mtx.gz"))))
  features_ref = read.csv(gzfile(file.path(comp_dir, paste0(organ, "_features.tsv.gz"))),
                          header = FALSE, stringsAsFactors = FALSE)[,1]
  features_ref = make.unique(features_ref, sep = "_dup")
  barcodes_ref = read.csv(gzfile(file.path(comp_dir, paste0(organ, "_barcodes.tsv.gz"))),
                          header = FALSE, stringsAsFactors = FALSE)[,1]
  metadata_ref = read.csv(file.path(comp_dir, paste0(organ, "_metadata.csv")), stringsAsFactors = FALSE)
  
  if (!"donor_id" %in% colnames(metadata_ref)) {
    metadata_ref$donor_id = "1"
  }
  
  rownames(counts_ref) = barcodes_ref
  colnames(counts_ref) = features_ref
  
  if(organ == "LI"){
    counts_ref = counts_ref[metadata_ref$cell,]
    row.names(metadata_ref) = metadata_ref$cell
  }
  
  counts_ST = Matrix::readMM(gzfile(file.path(organ_dir, organ, paste0(organ, "_X.mtx.gz"))))
  features_ST = read.csv(
    gzfile(file.path(organ_dir, organ, paste0(organ, "_features.tsv.gz"))),
    header = FALSE, stringsAsFactors = FALSE
  )[, 1]
  barcodes_ST = read.csv(
    gzfile(file.path(organ_dir, organ, paste0(organ, "_barcodes.tsv.gz"))),
    header = FALSE, stringsAsFactors = FALSE
  )[, 1]
  metadata_ST = read.csv(
    file.path(organ_dir, organ, paste0(organ, "_metadata.csv")),
    stringsAsFactors = FALSE
  )
  
  rownames(counts_ST) = barcodes_ST
  colnames(counts_ST) = features_ST
  
  # Create seurat objs (features x cells)
  
  counts_ref = Matrix::t(counts_ref)
  counts_ST  = Matrix::t(counts_ST)
  
  seurat_ref = CreateSeuratObject(counts_ref, meta.data = metadata_ref)
  seurat_ST_all  = CreateSeuratObject(counts_ST,  meta.data = metadata_ST)
  
  spots_keep = seurat_ST_all@meta.data %>%
    rownames_to_column("spot_names") %>%
    filter(Sample == sample) 
  
  seurat_ST = subset(seurat_ST_all, cells = spots_keep$spot_names)

  # Process seurat objects
  
  cells_keep = seurat_ref@meta.data %>%
    rownames_to_column("cell_names") %>%
    filter(!is.na(CellType_Deconvolution),
           !CellType_Deconvolution %in% c("unknown","Unknown","Other","other")) %>%
    group_by(CellType_Deconvolution) %>%
    group_modify(~ dplyr::slice_sample(.x, n = min(1000, nrow(.x)))) %>%
    mutate(cell_count = n()) %>%
    ungroup() %>%
    filter(cell_count>24)

  seurat_ref = subset(seurat_ref, cells = cells_keep$cell_names)
  
  min_cells = 5
  assay_ref = DefaultAssay(seurat_ref)
  assay_st  = DefaultAssay(seurat_ST)
  
  get_counts_mat = function(obj, assay){
    if (inherits(obj[[assay]], "Assay5")) {
      lyr = SeuratObject::Layers(obj[[assay]])
      lyr_counts = if ("counts" %in% lyr) "counts" else lyr[grepl("^counts", lyr)][1]
      if (is.na(lyr_counts)) stop("No counts-like layer found in Assay5")
      SeuratObject::LayerData(obj[[assay]], layer = lyr_counts)
    } else {
      GetAssayData(obj, assay = assay, slot = "counts")
    }
  }
  
  counts_ref = get_counts_mat(seurat_ref, assay_ref)
  counts_st  = get_counts_mat(seurat_ST,  assay_st)
  
  genes_ref = rownames(counts_ref)[Matrix::rowSums(counts_ref > 0) >= min_cells]
  genes_st  = rownames(counts_st)[Matrix::rowSums(counts_st  > 0) >= min_cells]
  
  common_genes = intersect(genes_ref, genes_st)
  
  seurat_ref = subset(seurat_ref, features = common_genes)
  seurat_ST = subset(seurat_ST, features = common_genes)
  
  
  coords_ST = data.frame(
    y = seurat_ST@meta.data$`y_scaled_image`,
    x = seurat_ST@meta.data$`x_scaled_image`
  )
  rownames(coords_ST) = row.names(seurat_ST@meta.data)
  



  # Build RCTD objects
  
  seurat_ref@meta.data$CellType_Deconvolution = factor(seurat_ref@meta.data$CellType_Deconvolution, levels = unique(seurat_ref@meta.data$CellType_Deconvolution))
  ref_counts = seurat_ref@assays$RNA$counts          # genes x cells
  ref_types  = seurat_ref@meta.data$CellType_Deconvolution
  names(ref_types) = colnames(ref_counts)
  ref_nUMI  = Matrix::colSums(ref_counts)
  
  sp_counts = seurat_ST@assays$RNA$counts             # genes x spots
  sp_nUMI   = Matrix::colSums(sp_counts)
  
  reference   = spacexr::Reference(ref_counts, ref_types, ref_nUMI)
  spatialRNA  = spacexr::SpatialRNA(counts = sp_counts, coords = coords_ST, nUMI = sp_nUMI)
  
  rctd = create.RCTD(
    spatialRNA = spatialRNA,
    reference = reference,
    max_cores = max(1, parallel::detectCores() - 1)
  )
  
  # Run RCTD
  
  rctd = spacexr::run.RCTD(rctd, doublet_mode = "doublet")
  
  # Process output
  out_dir = file.path(rctd_dir, organ)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  saveRDS(rctd, file = file.path(out_dir, paste0(organ, "_", sample, "_RCTD_object.rds")))
  

  results = rctd@results  # list with weights
  W = results$weights     # spots x celltypes
  
  organ_prop = as.data.frame(W) %>%
    tibble::rownames_to_column("barcode") %>%
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "CellType_Deconvolution",
      values_to = "proportion"
    )
  
  organ_prop = coords_ST %>%
    tibble::rownames_to_column("barcode") %>%
    left_join(organ_prop, by = "barcode") %>%
    mutate(Organ = organ)
  
  
  readr::write_csv(organ_prop, file.path(out_dir, paste0(organ, "_", sample, "_RCTD_proportions.csv")))
  
  invisible(list(RCTD = rctd, proportions = organ_prop))
}


# Run RCTD ------------------

samples = c("CTRL_1", "CTRL_2", "LPS_1", "LPS_2")

walk(samples, function(sample){
  Run_RCTD(organ, sample = sample)
})

