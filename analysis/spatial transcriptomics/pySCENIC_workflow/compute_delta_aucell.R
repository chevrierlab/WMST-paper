# SETUP --------------------------------------------------------------------------------------------------

source("E:/Coding/Scripts/R_Sourcing/R_Master_Source.R")
Combinations = read.csv2("E:/Coding/Scripts/ArraySeq/Projects/WMST_Revisions/PyScenic_Project/Run/IRF1_Stat1/Trimming_Plots/output/dfs/OrganxCore_Combinations_df_toKeep_afterCountThreshold_for_pyscenic.csv")
Metadata = read.csv2("Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/IHC/dfs/adata_obs.csv")


# COMPILE AUCell DATA -------------------------------------------------------------------


files = list.files(
  "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/Celltypes/Analysis/IRF1_Stat1_AuCell/regulon_AUCELL/1000",
  pattern = "\\.csv$",
  recursive = T,
  full.names = T
)

type_to_core = c(
  adipocytes       = "adipocytes",
  blood            = "blood & immune cell",
  endothelial_cell = "endothelial cell",
  epithelial_cell  = "epithelial cell",
  glial_cell       = "glial cell",
  mesenchymal_cell = "mesenchymal cell",
  muscle_cell      = "muscle cell",
  precursor_cell   = "precursor cell",
  stem_cell        = "stem cell"
)

summary_scores = map_dfr(files, function(path){
  
  auc_max = dirname(dirname(path)) %>% basename()
  base = basename(path)
  tf = basename(dirname(path))
  
  cell_unformatted = str_split(basename(path), "__")[[1]][2]
  organ =  str_split(basename(path), "__")[[1]][1]
  
  if(cell_unformatted %in% names(type_to_core)){
    cell = type_to_core[cell_unformatted]
    combo = paste(organ, cell, sep = "|")
    TF_combo = paste(tf, combo, sep = "|")
    
    if(combo %!in% Combinations$combo){return(NULL)}
  } else {return(NULL)}
  
  df = read.csv(path)
  
  if(nrow(df) != sum(df$spot %in% Metadata$spot)) {print("ERROR")}
  
  df = df %>%
    filter(spot %in% Metadata$spot) %>%
    left_join(
      .,
      Metadata %>% dplyr::rename(CellType_Core_Original = CellType_Core) %>% select(spot, Treatment, Sample, CellType_Core_Original, CellType_Broad),
      by = "spot"
    ) %>%
    mutate(
      TF = tf,
      auc_max = auc_max
    ) 

  CellType_Core = unique(df$CellType_Core)
  
  printg("loading [{auc_max}] {tf} | {organ} | {CellType_Core}  ~   {basename(path)}")
  
  df %>%
    group_by(auc_max_rank, TF, Organ, CellType_Core_Original, CellType_Core, Treatment, Sample) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(Type = "Core")
  
})

df_AUC_core = summary_scores %>%
  group_by(auc_max_rank, TF, Organ, CellType_Core_Original, CellType_Core, Treatment) %>%
  summarise(mean_val = mean(mean_score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Treatment,
    values_from = mean_val
  ) %>%
  mutate(delta_LPS = LPS - Control) %>% 
  as.df %>%
  mutate(combo = paste(Organ, CellType_Core_Original , sep = "|"),
         TF_combo = paste(TF, combo, sep = "|")) %>%
  select(-CellType_Core) %>%
  dplyr::rename(CellType_Core = CellType_Core_Original) %>%
  mutate(
    TF = factor(TF, levels = c("Irf1","Stat1", "Shared")),
    Organ = as.character(Organ),
    CellType_Core = as.character(CellType_Core)
  ) 

write.csv2(
  df_AUC_core,
  "E:/Coding/Scripts/ArraySeq/Projects/WMST_Revisions/PyScenic_Project/Clean_Run/IRF1_Stat1/Finalizing_Figures/output/AUC/AUC_Score_deltaLPS_aucMax_1000_STAT!_IRF1_Shared.csv",
  project = "WMST",
  topic = "WMST - Pyscenic",
  desc = "AUCell score summary for organ x core combinations. AUC max of 1000, for IRF1, Stat1, and shared regulons."
)
