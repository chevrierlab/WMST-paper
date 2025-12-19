# SETUP --------------------------------------------------

source("E:/Coding/Scripts/R_Sourcing/R_Master_Source.R")
source("E:/Coding/Scripts/ArraySeq/Projects/WMST_Revisions/PyScenic_Project/Run/IRF1_Stat1/parse_regulons_functions.R")

Combinations = read.csv("E:/Coding/Scripts/ArraySeq/Projects/WMST_Revisions/PyScenic_Project/Run/IRF1_Stat1/Trimming_Plots/output/dfs/OrganxCore_Combinations_df_toKeep_afterCountThreshold_for_pyscenic.csv")
root_dir = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/Celltypes/output/core"


# PARSE MOTIFS  -----------------------------------------

reg_files = list.files(root_dir, pattern = "\\.regulons\\.csv$", recursive = TRUE, full.names = TRUE) %>%
  get_reg_files(., c("Irf1", "Stat1"))


regulons_df = reg_files %>%
  map(parse_regulons_Irf1Stat1) %>%
  list_rbind() %>%
  dplyr::rename(gene = genes)

Stat1_regulons_df = regulons_df %>% filter(TF == "Stat1")
Irf1_regulons_df = regulons_df %>% filter(TF == "Irf1")
Shared_regulons_df = regulons_df %>%
  group_by(Organ, CellType_Core, gene) %>%
  filter(n()>1) %>%
  ungroup() %>%
  distinct(Organ, CellType_Core, gene, .keep_all = T) %>%
  mutate(TF = "Shared")

regulons_df = regulons_df %>%
  rbind(Shared_regulons_df)

regulons_df %>%
  distinct(TF, Organ, CellType_Core) %>%
  dplyr::count(TF)

Regulon_dir = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/Celltypes/Analysis/IRF1_Stat1_Regulons"

write.csv(regulons_df,
          file = file.path(Regulon_dir,"IRF1_Stat1_Regulons_OrganxCore.csv"))
write.csv(Stat1_regulons_df,
          file = file.path(Regulon_dir,"Stat1_Regulons_OrganxCore.csv"))
write.csv(Irf1_regulons_df,
          file = file.path(Regulon_dir,"IRF1_Regulons_OrganxCore.csv"))
write.csv(Shared_regulons_df,
          file = file.path(Regulon_dir,"Shared_Regulons_OrganxCore.csv"))

