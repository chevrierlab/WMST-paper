# SETUP -------------------

# Stat1/Irf1/shared regulon genes ∩ sig DEGs for core cell types within organs; cluster by membership pattern across core cell types; annotate clusters with GO BP + WikiPathways.

.libPaths(c("C:/Users/mhcle/Documents/R/4.3_upgradedRCPP",
            "C:/Users/mhcle/AppData/Local/R/win-library/4.3",
            .Library))
install.packages("shadowtext", repos = "https://cloud.r-project.org", type = "binary")
library(clusterProfiler)


source("E:/Coding/Scripts/R_Sourcing/R_Master_Source.R")
source("E:/Coding/Scripts/ArraySeq/Projects/WMST_Revisions/PyScenic_Project/Clean_Run/IRF1_Stat1/Finalizing_Figures/Finalize_Organ_HMs/Functions/Organ_HM_Functions.R")

DEG_dir = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/scenic/Celltypes/Analysis/tracking/DEGs_associated_withRegulons/Rerunning_DEGs/Final_DEG_Lists"
DEGs_Organ_Core = read.csv2(file.path(DEG_dir,"DEGs_Organ_Core.csv"))
regulons_df = read.csv2("E:/Coding/Scripts/ArraySeq/Projects/WMST_Revisions/PyScenic_Project/Clean_Run/IRF1_Stat1/Trimming_Plots/output/dfs/Regulons_df_IRF1_Stat1_afterCountThreshold.csv")

# HEATMAPPER ---------------------------------------------------------------------

organs = unique(regulons_df$Organ)

walk(organs, function(organ){
  
  organ_full = Metadata %>%
    filter(Organ == organ) %>%
    pull(Organ_Full_Name) %>%
    unique
  TAG = organ
  
  # Filter dfs
  regulon_f1 = regulons_df %>%
    filter(Organ == organ) %>%
    mutate(Organ_x_Core = paste(Organ, CellType_Core, sep = "|"))
  
  Organ_Core_DEG_f = DEGs_Organ_Core %>%
    filter(Organ_x_Core %in% (regulon_f1$Organ_x_Core %>% unique())) %>%
    filter(gene %in% (regulon_f1$gene %>% unique()))%>%
    filter(pval_adj < 0.05, abs(logfc) > 0.25)
  
  regulon_f = Organ_Core_DEG_f %>%
    select(gene, Organ_x_Core) %>%
    left_join(.,
              regulon_f1,
              by = c("gene", "Organ_x_Core"))
  
  Regulon_genes = length(unique(regulon_f$gene))
  Organ_Core_genes = length(unique(Organ_Core_DEG_f$gene))
  
  printg("[{TAG}] {Regulon_genes} gene across regulons. 
         {percenter(Organ_Core_genes, Regulon_genes)}% overlap with Core DEGs.
         Expected is 100%")
  
  
  set.seed(101)
  
  
  # Construct heatmaps ----------------------------------------------------------------------------------
  
  # Set save dir
  base_dir = "E:/Coding/Scripts/ArraySeq/Projects/WMST_Revisions/PyScenic_Project/Clean_Run/IRF1_Stat1/Finalizing_Figures/output/Heatmap_Organ/"
  out_dir_bucket =  folder.path(base_dir, organ, "Bucket")
  out_dir_pathways_bucket = folder.path(out_dir_bucket, "pathway_dfs")
  
  HM_Membership_bucket = Membership_HM_buckets(
    reg_df=regulon_f,
    core_df = Organ_Core_DEG_f
  )
  
  HM_CoreDE = Core_DEG_HM(
    core_df=Organ_Core_DEG_f, 
    reg_df=regulon_f, 
    row_order = HM_Membership_bucket$rows, 
    col_order = HM_Membership_bucket$cols 
  )
  
  HM_CoreLFC = Core_LFC_HM(
    core_df=Organ_Core_DEG_f, 
    reg_df=regulon_f, 
    row_order = HM_Membership_bucket$rows, 
    col_order = HM_Membership_bucket$cols 
  )
  
  
  
  # Pathway Enrichment
  
  DE_all = DEGs_Organ_Core %>%
    filter(Organ_x_Core %in% (regulons_df$combo %>% unique())) %>%
    filter(pval_adj < 0.05, abs(logfc) > 0.25)
  
  
  regulon_all = inner_join(
    DE_all %>% select(gene, Organ_x_Core),
    regulons_df %>% dplyr::rename(Organ_x_Core = combo),
    by = c("gene", "Organ_x_Core")
  )
  
  universe_genes = unique(regulon_all$gene) 
  
  
  pathway_df = map_dfr(unique(HM_Membership_bucket$cluster_df$cluster), function(clust){
    genes = HM_Membership_bucket$cluster_df %>% filter(cluster == clust) %>% pull(gene) %>% as.vector
    
    pathways = pathway_enrich_UNIVERSE(genes = genes, cluster = clust, universe = universe_genes)
    
    write.csv(pathways$df,
              file = file.path(out_dir_pathways_bucket, paste0("Cluster_", clust,"_Pathways.csv")),
              row.names = F)
    
    pathways$more_selected %>%
      mutate(
        cluster = clust
      )
  })
  
  
  # Annotate heatmaps
  
  HM_Membership_ann = Membership_HM_anno(
    Cluster_df = HM_Membership_bucket$cluster_df,
    matrix = HM_Membership_bucket$matrix, 
    col_order = HM_Membership_bucket$cols,
    pathway_df = pathway_df
  )
  
  
  Save_Final_Heatmaps(
    HM1 = HM_Membership_ann, 
    HM2 = HM_CoreDE, 
    HM3 = HM_CoreLFC, 
    title = paste0(organ_full, ": Membership-bucketed clustering, Wiki and GO BP top 10"), 
    out_dir = out_dir_bucket, 
    h = 10, 
    w = 10
  )
  
  
  HM_Membership_ann_TopGo = Membership_HM_anno(
    Cluster_df = HM_Membership_bucket$cluster_df,
    matrix = HM_Membership_bucket$matrix, 
    col_order = HM_Membership_bucket$cols,
    pathway_df = pathway_df,
    type_pathway = "GO_BP",
    pathway_count = 3
  )
  
  HM_Membership_ann_Top_Wiki = Membership_HM_anno(
    Cluster_df = HM_Membership_bucket$cluster_df,
    matrix = HM_Membership_bucket$matrix, 
    col_order =  HM_Membership_bucket$cols,
    pathway_df = pathway_df,
    type_pathway = "Wiki",
    pathway_count = 3
  )
  
  Save_Final_Heatmaps(
    HM1 = HM_Membership_ann_TopGo, 
    HM2 = HM_CoreDE, 
    HM3 = HM_CoreLFC, 
    title = paste0(organ_full, ": Membership-bucketed clustering, GO BP top 3"), 
    out_dir = out_dir_bucket, 
    h = 10, 
    w = 10
  )
  
  Save_Final_Heatmaps(
    HM1 = HM_Membership_ann_Top_Wiki, 
    HM2 = HM_CoreDE, 
    HM3 = HM_CoreLFC, 
    title = paste0(organ_full, ": Membership-bucketed clustering, Wiki top 3"), 
    out_dir = out_dir_bucket, 
    h = 10, 
    w = 10
  )
  
  
  
}