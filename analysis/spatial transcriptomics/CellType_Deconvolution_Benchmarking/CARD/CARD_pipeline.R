# SETUP -----------------------------------------------------

user_lib = file.path(Sys.getenv("HOME"), "R", paste0(R.version$major, ".", sub("\\..*$","", R.version$minor)))
library(TOAST, lib.loc = user_lib)
library(MuSiC, lib.loc = user_lib)
library(sf); library(V8);
library(CARD)
library(RhpcBLASctl)

source("/project/nchevrier/projects/clevenger1/Projects/Resources/Functions_Sourcing/R_Master_Source.R")

organs = c("SP", "LU", "LI", "KI")
organ  = organs[as.integer(Sys.getenv("ARRAY_ID")) + 1]

print(organ)

# SET THREADS --------------------------------

blas_set_num_threads(8)
omp_set_num_threads(8)
Sys.setenv(
  OMP_NUM_THREADS        = "8",
  OPENBLAS_NUM_THREADS   = "8",
  MKL_NUM_THREADS        = "8",
  VECLIB_MAXIMUM_THREADS = "8",
  NUMEXPR_NUM_THREADS    = "8"
)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(8)


# CARD Pipeline -----------------------------------------------------------------------------------


Run_CARD = function(
    organ,
    organ_dir = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Organ_exports" %>% spath,
    base_dir  = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Reference_data" %>% spath,
    card_dir = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/CARD_output" %>% spath,
    sample
){
  
  # Load components
  print(paste(organ, ": loading reference"))
  
  comp_dir = file.path(base_dir, organ, "components")
  counts_ref = Matrix::readMM(gzfile(file.path(comp_dir, paste0(organ, "_X_counts.mtx.gz"))))
  features_ref = read.csv(gzfile(file.path(comp_dir, paste0(organ, "_features.tsv.gz"))),
                          header = FALSE, stringsAsFactors = FALSE)[,1]
  barcodes_ref = read.csv(gzfile(file.path(comp_dir, paste0(organ, "_barcodes.tsv.gz"))),
                          header = FALSE, stringsAsFactors = FALSE)[,1]
  metadata_ref = read.csv(file.path(comp_dir, paste0(organ, "_metadata.csv")), stringsAsFactors = FALSE)
  
  rownames(counts_ref) = barcodes_ref
  colnames(counts_ref) = features_ref
  
  if (organ == "LI") {
    counts_ref = counts_ref[metadata_ref$cell, ]
    row.names(metadata_ref) = metadata_ref$cell
  }
  
  print(paste(organ, ": loading ST"))
  
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
  
  # Create seurat objs (Seurat expects features x cells)
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
  
  # Create Card object
  coords_ST = data.frame(
    y = seurat_ST@meta.data$`y_scaled_image`,
    x = seurat_ST@meta.data$`x_scaled_image`
  )
  rownames(coords_ST) = row.names(seurat_ST@meta.data)
  
  seurat_ref@meta.data$CellType_Deconvolution = factor(seurat_ref@meta.data$CellType_Deconvolution, levels = unique(seurat_ref@meta.data$CellType_Deconvolution))
  if (!"donor_id" %in% colnames(seurat_ref@meta.data)) {
    seurat_ref@meta.data$donor_id = "1"
  }
  seurat_ref@meta.data$donor_id = factor(seurat_ref@meta.data$donor_id, levels = unique(seurat_ref@meta.data$donor_id))
  
  CARD_obj = createCARDObject(
    sc_count = seurat_ref@assays$RNA$counts,
    sc_meta  = seurat_ref@meta.data,
    spatial_count     = seurat_ST@assays$RNA$counts,
    spatial_location  = coords_ST,
    ct.varname        = "CellType_Deconvolution",
    ct.select         = unique(seurat_ref@meta.data$CellType_Deconvolution),
    sample.varname    = "donor_id",
    minCountGene = 5,
    minCountSpot = 5
  )
  
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  
  out_dir = file.path(card_dir, organ)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  saveRDS(CARD_obj, file = file.path(out_dir, paste0(organ, "_", sample, "_CARD_object.rds")))
  
  organ_prop = as.data.frame(CARD_obj@Proportion_CARD) %>%
    tibble::rownames_to_column("barcode") %>%
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "CellType_Deconvolution",
      values_to = "proportion"
    )
  
  spatial_locs = as.data.frame(CARD_obj@spatial_location) %>%
    tibble::rownames_to_column("barcode") %>%
    left_join(organ_prop, by = "barcode") %>%
    mutate(Organ = organ)
  
  readr::write_csv(spatial_locs, file.path(out_dir, paste0(organ, "_", sample, "_CARD_proportions.csv")))
  
  invisible(list(CARD_obj = CARD_obj, proportions = spatial_locs))
}


# Run CARD -----------------------------------------------------------------------------------

samples = c("CTRL_1", "CTRL_2", "LPS_1", "LPS_2")

walk(samples, function(sample){
  Run_CARD(organ, sample = sample)
})

