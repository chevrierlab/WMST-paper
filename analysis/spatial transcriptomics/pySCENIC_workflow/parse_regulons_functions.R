# SETUP --------------------------------------------------------------------------


cc_to_core = c(
  "adipocytes"        = "adipocytes",
  "blood_immune_cell" = "blood & immune cell",
  "endothelial_cell"  = "endothelial cell",
  "epithelial_cell"   = "epithelial cell",
  "glandular_cell"    = "glandular cell",
  "glial_cell"        = "glial cell",
  "mesenchymal_cell"  = "mesenchymal cell",
  "muscle_cell"       = "muscle cell",
  "neuronal_cell"     = "neuronal cell",
  "precursor_cell"    = "precursor cell",
  "stem_cell"         = "stem cell"
)

# FUNCTIONS ----------------------------------------------------------------------



get_reg_files = function(reg_files, TF){
  reg_files[
    vapply(reg_files, function(f) {
      df = tryCatch(
        read.csv(f, header = TRUE, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      !is.null(df) &&
        nrow(df) > 1 &&
        any(df[[1]] %in% TF)
    }, logical(1))
  ] 
}



filter_motif = function(df, tf, Organ, CellType_Core){
  
  TF_df =  df %>% filter(TF == tf)
  motif_count = TF_df %>% nrow
  
  df_filtered = TF_df %>%
    filter(NES >= 3.2, Motif_Similarity_Qvalue <=  0.001)
  
  if(nrow(df_filtered)>0){
    TF_motifs = map(seq_len(nrow(df_filtered)), function(i){
      gene_list = df_filtered[["Genes"]][i]
      str_split(gene_list, ",\\s*") %>% unlist() %>% unique()
    })
    
    TF_regulon = tibble(set_id = seq_along(TF_motifs), genes = TF_motifs) %>%
      unnest(genes) %>%
      distinct(set_id, genes) %>%
      dplyr::count(genes, name = "n_sets") %>%
      arrange(desc(n_sets)) %>%
      as.df %>%
      filter(n_sets >= 2) %>%
      pull(genes)
    
    TF_df_final =  data.frame(
      TF = tf,
      Organ = Organ,
      CellType_Core = CellType_Core,
      genes = TF_regulon,
      n_motifs = length(TF_motifs),
      n_genes = TF_regulon %>% length()
    )
  
     return(TF_df_final)
  
    } else{
    
      return(NULL)

  }
  
  
}


parse_regulons_Irf1Stat1 = function(f, cc_to_core = cc_to_core, Combinations = Combinations) {
  
  
  filename = basename(dirname(f))
  Organ = str_split(filename, "_")[[1]][1]
  CellType_Core_unformatted = gsub(
    paste0(Organ, "_"), "", filename
  )
  CellType_Core = cc_to_core[CellType_Core_unformatted] %>% as.character()
  
  combo = paste(Organ, CellType_Core, sep = "|")
  
  if(combo %!in% Combinations$combo){
    return(NULL)
  } 
  
  regulon_df = f %>% read.csv %>% as.df
  row_num = which(regulon_df[[1]] %in% c("Irf1", "Stat1"))[1]
  regulon_df = regulon_df[row_num:nrow(regulon_df),]
  colnames(regulon_df) = c("TF", "MotifID", "AUC", "NES", "Motif_Similarity_Qvalue",	"Orthologous_Identity",	"Annotation",
                           "Context",	"TargetGenes",	"RankAtMax") 
  
  
  
  regulon_df_processed = regulon_df %>%
    mutate(Genes = map_chr(TargetGenes, function(x) {
      if (is.na(x) || !nzchar(x)) return("")
      m = str_match_all(x, "\\('([^']+)'\\s*,")[[1]]  # capture the gene name before the comma
      if (nrow(m) == 0) return("")
      paste(m[,2], collapse = ", ")
    })) %>%
    mutate(
      NES = as.numeric(NES),
      Motif_Similarity_Qvalue = as.numeric(Motif_Similarity_Qvalue)
    )
  
  
  df_summary = rbind(
    filter_motif(
      df = regulon_df_processed,
      tf = "Stat1",
      Organ = Organ,
      CellType_Core = CellType_Core
    ),
    filter_motif(
      df = regulon_df_processed,
      tf = "Irf1",
      Organ = Organ,
      CellType_Core = CellType_Core
    )
  )
  
  return(df_summary)
  
}




